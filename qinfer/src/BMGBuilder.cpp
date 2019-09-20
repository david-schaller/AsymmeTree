#include <algorithm>
#include <iterator>
#include <iostream>
#include <random>

#include "BMGBuilder.h"

/******************************************************************************
                               common functions
******************************************************************************/

void
BMGBuilder::buildBMG(){
  // add all gene pointers to the graph
  for(Gene& gene : m_ptrS->getGenes()){
    m_bmg.addNode(&gene);
  }

  // use the simple epsilon method if quartets are disabled
  if(m_ptrParam->quartetsDisabled()){
    if(m_benchmark) m_benchmark->startEpsilon();
    epsilonMethod();
    if(m_benchmark) m_benchmark->endEpsilon();
    return;
  }
  if(m_benchmark) m_benchmark->startEpsilon();

  // build a matrix to check whether y is a best match candidate for x
  if(m_ptrParam->restrictedY()){
    if(m_benchmark) m_benchmark->startEpsilon();
    buildCandidateMatrix();
    if(m_benchmark) m_benchmark->endEpsilon();
  }

  // main algorithms
  if(!m_ptrParam->relativeOutgroups()){
    if(m_benchmark) m_benchmark->startBuildBMG();
    buildRootOutgroups();
    if(m_benchmark) m_benchmark->endBuildBMG();
  } else {
    m_outgroupChoice.initialize();
    if(m_benchmark) m_benchmark->startBuildBMG();
    buildRelativeOutgroups();
    if(m_benchmark) m_benchmark->endBuildBMG();
  }

}

void
BMGBuilder::printBMG(){
  for(const auto& vNbrsPair : m_bmg.getAdjacency()){
    if(vNbrsPair.second.size() > 0){
      for(const Gene* nbr : vNbrsPair.second){
        std::cout << vNbrsPair.first->getIdentifier()
                  << " "
                  << nbr->getIdentifier()
                  << std::endl;
      }
    }
  }
}

void
BMGBuilder::buildCandidateMatrix(){
  size_t dim {m_ptrS->getGenes().size()};
  m_bmCandidates.initMatrix(dim, 0);
  double threshold = 1.0 + m_ptrParam->getEpsilon();

  for(size_t i = 0; i < dim; ++i){
    auto speciesMap = std::unordered_map<std::string, double>();
    const auto& species_i = m_ptrS->getGeneSpecies(i);

    for(size_t j = 0; j < dim; ++j){
      const auto& species_j = m_ptrS->getGeneSpecies(j);
      if(species_i != species_j){
        if(speciesMap.find(species_j) == speciesMap.end()){
          // initialize minimum if key did not exist so far
          speciesMap[species_j] = m_ptrS->getDistance(i, j);
        } else {
          // update the minimum
          if(m_ptrS->getDistance(i, j) < speciesMap.at(species_j)){
            speciesMap.at(species_j) = m_ptrS->getDistance(i, j);
          }
        }
      }
    }

    for(size_t j = 0; j < dim; ++j){
      const auto& species_j = m_ptrS->getGeneSpecies(j);
      if( (species_i != species_j) &&
          (m_ptrS->getDistance(i, j) <= threshold * speciesMap.at(species_j))){
        m_bmCandidates.at(i,j) = 1;
      }
    }
  }

  // // test output of the matrix
  // for(size_t i = 0; i < dim; ++i){
  //   for(size_t j = 0; j < dim; ++j){
  //       std::cout << m_bmCandidates.at(i,j) << "  ";
  //   }
  //   std::cout << std::endl;
  // }
}

/******************************************************************************
                          Extended Best Hits method
******************************************************************************/
void
BMGBuilder::epsilonMethod(){
  size_t dim {m_ptrS->getGenes().size()};
  double threshold = 1.0 + m_ptrParam->getEpsilon();

  for(size_t i = 0; i < dim; ++i){
    auto speciesMap = std::unordered_map<std::string, double>();
    const auto& species_i = m_ptrS->getGeneSpecies(i);

    for(size_t j = 0; j < dim; ++j){
      const auto& species_j = m_ptrS->getGeneSpecies(j);
      if(species_i != species_j){
        if(speciesMap.find(species_j) == speciesMap.end()){
          // initialize minimum if key did not exist so far
          speciesMap[species_j] = m_ptrS->getDistance(i, j);
        } else {
          // update the minimum
          if(m_ptrS->getDistance(i, j) < speciesMap.at(species_j)){
            speciesMap.at(species_j) = m_ptrS->getDistance(i, j);
          }
        }
      }
    }

    for(size_t j = 0; j < dim; ++j){
      const auto& species_j = m_ptrS->getGeneSpecies(j);
      if( (species_i != species_j) &&
          (m_ptrS->getDistance(i, j) <= threshold * speciesMap.at(species_j))){
        m_bmg.addEdge(m_ptrS->getGenePtr(i), m_ptrS->getGenePtr(j));
      }
    }
  }
}

/******************************************************************************
                              OUTGROUG METHOD I
                   (choose outgroups w.r.t the root of S)
******************************************************************************/

std::vector<Gene*>
BMGBuilder::chooseOutgroups(const std::vector<Gene*>& outgroupCandidates){
  auto outgroups = std::vector<Gene*>();

  size_t limit = std::min(m_ptrParam->getOutgroupLimit(), outgroupCandidates.size());
  std::sample(outgroupCandidates.begin(), outgroupCandidates.end(),
              std::back_inserter(outgroups),
              limit,
              std::mt19937{std::random_device{}()});

  // for (auto genePtr : outgroups)
  //     std::cout << "outgroup" << genePtr->getIdentifier() << std::endl;

  return outgroups;
}

void
BMGBuilder::buildRootOutgroups(){

  for(Gene& x : m_ptrS->getGenes()){
    const std::vector<Gene*>& outgroupCandidates = m_ptrS->getOutgroups(&x);

    if(m_benchmark) m_benchmark->startChooseOutgroups();
    std::vector<Gene*> outgroupsZ = chooseOutgroups(outgroupCandidates);
    if(m_benchmark) m_benchmark->endChooseOutgroups();

    for(const std::string& speciesY : m_ptrS->getSpeciesSubtree(x.getSubtree())){

      // skip the own species
      if(x.getSpecies() == speciesY){
        continue;
      }

      // build the gene set Y
      const std::vector<Gene*>& allGenesY = m_ptrS->getSpeciesGenes(speciesY);
      std::vector<Gene*> genesY;
      if(!m_ptrParam->restrictedY()){
        genesY = std::vector<Gene*>(allGenesY);
      } else {
        genesY = std::vector<Gene*>();
        for(Gene* genePtr : allGenesY){
          if(m_bmCandidates.at(x.getIndex(), genePtr->getIndex())){
            genesY.push_back(genePtr);
          }
        }
      }

      if(genesY.size() == 0){
        std::cerr << "Warning: Species without genes: "
                  << speciesY
                  << std::endl;

      } else if (genesY.size() == 1){
        // trivial case: only one gene in the species of Y
        m_bmg.addEdge(&x, genesY[0]);

      } else {
        // find best matches in gene set Y
        for(Gene* bm : m_quartets.findBestMatches(&x, genesY, outgroupsZ)){
          m_bmg.addEdge(&x, bm);
        }
      }
    }
  }
}

/******************************************************************************
                              OUTGROUG METHOD II
          (relative outgroups corrected with incongruent quartets)
******************************************************************************/
void
BMGBuilder::buildRelativeOutgroups(){

  for(Gene& x : m_ptrS->getGenes()){
    // const std::vector<Gene*>& outgroupCandidates = m_ptrS->getOutgroups(&x);
    // std::vector<Gene*> outgroupsZ = chooseOutgroups(outgroupCandidates);

    for(const std::string& speciesY : m_ptrS->getSpeciesSubtree(x.getSubtree())){

      // skip the own species
      if(x.getSpecies() == speciesY){
        continue;
      }

      // build the gene set Y
      const std::vector<Gene*>& allGenesY = m_ptrS->getSpeciesGenes(speciesY);
      std::vector<Gene*> genesY;
      if(!m_ptrParam->restrictedY()){
        genesY = std::vector<Gene*>(allGenesY);
      } else {
        genesY = std::vector<Gene*>();
        for(Gene* genePtr : allGenesY){
          if(m_bmCandidates.at(x.getIndex(), genePtr->getIndex())){
            genesY.push_back(genePtr);
          }
        }
      }

      if(genesY.size() == 0){
        std::cerr << "Warning: Species without genes: "
                  << speciesY
                  << std::endl;

      } else if (genesY.size() == 1){
        // trivial case: only one gene in the species of Y
        m_bmg.addEdge(&x, genesY[0]);

      } else {
        // find best matches in gene set Y
        if(m_benchmark) m_benchmark->startChooseOutgroups();
        std::vector<Gene*> outgroupsZ = m_outgroupChoice.getClosest(&x, genesY);
        if(m_benchmark) m_benchmark->endChooseOutgroups();
        for(Gene* bm : m_quartets.findBestMatches(&x, genesY, outgroupsZ)){
          m_bmg.addEdge(&x, bm);
        }
      }
    }
  }
}
