#include <fstream>
#include <iostream>
#include <set>
#include <regex>

#include "Scenario.h"
#include "Tree.h"

void
Scenario::addGenes(std::deque<Gene>& g)
{
  m_genes.clear();

  m_genes.insert(m_genes.begin(), std::make_move_iterator(g.begin()),
                 std::make_move_iterator(g.end()));
  rebuildMap();
}

const std::vector<Gene*>&
Scenario::getOutgroups(Gene* genePtr) const {
  return m_outgroups[genePtr->getSubtree()];
}

const std::vector<std::string>&
Scenario::getSpeciesSubtree(size_t subtreeIdx) const {
  return m_STreeSubtrees[subtreeIdx];
}

std::vector<std::string>
Scenario::split(std::string s){
  std::regex ws_re("\\s+");   // split at whitespace(s)
  std::vector<std::string> result{
      std::sregex_token_iterator(s.begin(), s.end(), ws_re, -1), {}
  };
  return result;
}

void
Scenario::rebuildMap()
{
  m_geneAssignments.clear();

  for(auto it = m_genes.begin(); it != m_genes.end(); ++it) {
    m_geneAssignments.insert(
      std::pair<std::string, Gene*>(it->getIdentifier(), &(*it))
    );
  }
}

void
Scenario::parseFiles(){
  parseDistanceMatrix(m_ptrParam->getMatrixFilename());
  parseSpeciesGenes(m_ptrParam->getSpeciesFilename());
  if(!m_ptrParam->quartetsDisabled()){
    parseSTreeSubtrees(m_ptrParam->getTreeFilename(), m_ptrParam->subtreeFiles());
  }
}

void
Scenario::parseDistanceMatrix(std::string filepath){
  std::ifstream filestream(filepath, std::ios::binary);

  if(!filestream) {
    throw std::runtime_error("Failed to open file " + filepath);
  }

  size_t lineCounter = 0;
  std::string line;
  while(std::getline(filestream, line)) {
    // read the first line which contains the dimension
    if(lineCounter == 0) {
      m_distanceMatrix.initMatrix(std::stoi(line), 0.0);
    }
    else {
      parseDistanceMatrixRow(line, lineCounter-1);
    }

    ++lineCounter;
  }

  rebuildMap();

  // // test output of the matrix
  // for(size_t i = 0; i < m_distanceMatrix.getDim(); ++i){
  //   for(size_t j = 0; j < m_distanceMatrix.getDim(); ++j){
  //       std::cout << m_distanceMatrix.at(i,j) << "  ";
  //   }
  //   std::cout << std::endl;
  // }
}

void
Scenario::parseDistanceMatrixRow(std::string row, size_t rowIdx){

  std::vector<std::string> items = split(row);
  size_t columnIdx = 0;

  for(auto& item : items){
    if(columnIdx == 0){
      m_genes.push_back(Gene(item, m_genes.size()));
    } else {
      m_distanceMatrix.at(rowIdx, columnIdx - 1) = std::stod(item);
    }
    ++columnIdx;
  }
}

void
Scenario::parseSpeciesGenes(std::string filepath){
  std::ifstream filestream(filepath, std::ios::binary);

  if(!filestream) {
    throw std::runtime_error("Failed to open file " + filepath);
  }

  auto v = std::vector<std::vector<double>>();
  auto g = std::deque<Gene>();
  std::string line;
  while(std::getline(filestream, line)) {
    parseSpeciesGenesLine(line);
  }
}

void
Scenario::parseSpeciesGenesLine(std::string line){
  auto g = std::vector<Gene*>();
  std::string species;

  std::vector<std::string> items = split(line);
  bool firstElement = true;

  for(auto& item : items){
    // the first element is the species identifier
    if(firstElement) {
      species = item;
      firstElement = false;
    } else {
      try {
        g.push_back(m_geneAssignments.at(item));
        m_geneAssignments.at(item)->setSpecies(species);
      } catch(const std::out_of_range&) {
        std::cerr << "WARNING: Species-to-genes-file contains additional entries: "
                  << item
                  << std::endl;
      }
    }
  }

  m_speciesGenes.insert(std::pair<std::string, std::vector<Gene*>>(species, g));
}

void
Scenario::parseSTreeSubtrees(std::string filepath, bool subtreeFiles){
  std::ifstream filestream(filepath, std::ios::binary);

  if(!filestream) {
    throw std::runtime_error("Failed to open file " + filepath);
  }

  std::string line;

  if(!subtreeFiles){
    std::getline(filestream, line);
    parseNewickAndCheck(line);
  } else {
    while(std::getline(filestream, line)) {
      parseSTreeSubtreeLine(line);
    }
  }

  // for(auto it = m_STreeSubtrees.begin(); it != m_STreeSubtrees.end(); ++it){
  //   for(auto it2 = (*it).begin(); it2 != (*it).end(); ++it2){
  //     std::cout << *it2 << "  ";
  //   }
  //   std::cout << std::endl;
  // }

  // check if all species with genes are in the subtree lists
  checkSpeciesAvailability();
  // build the corresponding lists of outgroup genes
  if(m_benchmark) m_benchmark->startOutgroupInit();
  buildOutgroupInfo();
  if(m_benchmark) m_benchmark->endOutgroupInit();
}

void
Scenario::parseNewickAndCheck(std::string line){
  m_STree = Tree::parseNewick(line);
  std::vector<std::vector<std::string>> subtrees = m_STree.speciesInSubtrees();
  for(auto& subtree : subtrees){
    m_STreeSubtrees.push_back(std::vector<std::string>());
    for(auto& species : subtree){
      if(m_speciesGenes.find(species) == m_speciesGenes.end()){
        std::cerr << "WARNING: Newick tree contains species without genes: "
                  << species
                  << " (omitted)"
                  << std::endl;
      } else {
        m_STreeSubtrees.back().push_back(species);
      }
    }
  }
}

void
Scenario::parseSTreeSubtreeLine(std::string line){
  m_STreeSubtrees.push_back(std::vector<std::string>());

  std::vector<std::string> items = split(line);

  for(auto& item : items){
    if(m_speciesGenes.find(item) == m_speciesGenes.end()){
      std::cerr << "WARNING: Species-subtrees-file contains species without genes: "
                << item
                << " (omitted)"
                << std::endl;
    } else {
      m_STreeSubtrees.back().push_back(item);
    }
  }
}

void
Scenario::checkSpeciesAvailability(){
  auto necessarySpecies = std::set<std::string>();
  for(auto& species : m_speciesGenes){
    necessarySpecies.insert(species.first);
  }

  for(auto& subtree : m_STreeSubtrees){
    for(auto& species : subtree){
      if(necessarySpecies.find(species) == necessarySpecies.end()){
        throw std::runtime_error("Could not assign species to subtree: " + species);
      } else {
        necessarySpecies.erase(species);
      }
    }
  }

}

void
Scenario::buildOutgroupInfo(){
  for(size_t i = 0; i < m_STreeSubtrees.size(); ++i){
    m_outgroups.push_back(std::vector<Gene*>());
    for(size_t j = 0; j < m_STreeSubtrees.size(); ++j){
      if(i != j){
        for(auto species : m_STreeSubtrees[j]){
          m_outgroups.back().insert(m_outgroups.back().end(),
                                    m_speciesGenes[species].begin(),
                                    m_speciesGenes[species].end());
        }
      }
    }

    for(const auto& species : m_STreeSubtrees[i]){
      for(auto genePtr : m_speciesGenes[species]){
        genePtr->setSubtree(i);
      }
    }
  }

  // for(auto it = m_outgroups.begin(); it != m_outgroups.end(); ++it){
  //   for(auto it2 = (*it).begin(); it2 != (*it).end(); ++it2){
  //     std::cout << (*it2)->getIdentifier() << "  ";
  //   }
  //   std::cout << std::endl;
  // }
}
