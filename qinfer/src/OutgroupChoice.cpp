#include <algorithm>
#include <random>

#include <iostream>

#include "OutgroupChoice.h"

void
OutgroupChoice::quicksortRow(size_t row, int l, int r){
  if(r <= l){
    return;
  }
  int i = l;
  int j = r;
  double piv = m_ptrS->getDistance(row, m_I.at(row, (l+r)/2));

  while(i <= j){

    while(m_ptrS->getDistance(row, m_I.at(row, i)) < piv){
      ++i;
    }
    while(m_ptrS->getDistance(row, m_I.at(row, j)) > piv){
      --j;
    }

    if(i <= j){
      size_t temp = m_I.at(row, i);
      m_I.at(row, i) = m_I.at(row, j);
      m_I.at(row, j) = temp;
      ++i;
      --j;
    }
  }

  quicksortRow(row, l, j);
  quicksortRow(row, i, r);
}

void
OutgroupChoice::buildIMatrix(){
  size_t dim {m_ptrS->getGenes().size()};
  m_I.initMatrix(dim, 0);

  for(size_t j = 1; j < dim; ++j){
    for(size_t i = 0; i < dim; ++i){
      m_I.at(i,j) = j;
    }
  }

  for(size_t row = 0; row < dim; ++row){
    quicksortRow(row, 0, dim-1);
  }

  // // test output of the matrix
  // for(size_t i = 0; i < dim; ++i){
  //   for(size_t j = 0; j < dim; ++j){
  //       std::cout << m_I.at(i,j) << "  ";
  //   }
  //   std::cout << std::endl;
  // }
}

void
OutgroupChoice::computeLcaS(){
  m_ptrS->getSTree().supplyLeaves();
  size_t nodeN = m_ptrS->getSTree().getNodeNumber();
  size_t leafN = m_ptrS->getSTree().getLeafNumber();

  m_lcaS.initMatrix(leafN, 0);
  m_subtreeSpecies = std::vector<std::vector<std::shared_ptr<TreeNode>>>(nodeN, std::vector<std::shared_ptr<TreeNode>>());
  m_subtreeGenes = std::vector<std::vector<Gene*>>(nodeN, std::vector<Gene*>());

  for(auto& v : m_ptrS->getSTree().getPreorder()){

    size_t vid = v->getNodeIdx();
    // m_nodeIdxToNode[vid] = v->getValue();                           // help map
    m_subtreeSpecies[vid].insert(m_subtreeSpecies[vid].end(),
                                 v->getLeaves().begin(),
                                 v->getLeaves().end());
    // std::cout << "N" << v->getValue();
    // for(auto species : m_subtreeSpecies[vid]){
    //   std::cout << "  " << species->getValue();
    // }
    // std::cout << std::endl;

    for(auto& leaf : m_subtreeSpecies[vid]){
      auto genes = m_ptrS->getSpeciesGenes(leaf->getValue());
      m_subtreeGenes[vid].insert(m_subtreeGenes[vid].end(),
                                 genes.begin(), genes.end());
    }

    auto& children = v->getChildren();
    if(!v->hasChildren()){
      m_lcaS.at(v->getLeafIdx(), v->getLeafIdx()) = vid;
      // fill the species identifier --> leaf index map
      m_speciesLeafIdx[v->getValue()] = v->getLeafIdx();

    } else if (children.size() >= 2){
      // for all combinations of 2 children of v
      for(size_t i = 0; i < children.size()-1; ++i){
        for(size_t j = i; j < children.size(); ++j){

          // for all combinations of 2 genes from the respective species
          for(auto& l1 : children[i]->getLeaves()){
            for(auto& l2 : children[j]->getLeaves()){
              m_lcaS.at(l1->getLeafIdx(), l2->getLeafIdx()) = vid;
              m_lcaS.at(l2->getLeafIdx(), l1->getLeafIdx()) = vid;
            }
          }
        }
      }
    }
  }
}

void
OutgroupChoice::computeOutgroups(){

  for(size_t i = 0; i < m_subtreeSpecies.size(); ++i){
    m_lcaOutgroups[i] = std::unordered_set<std::string>();
  }

  // the outgroups on the root of S are fixed
  bool rootDone = false;

  for(auto& v : m_ptrS->getSTree().getPreorder()){
    size_t vid = v->getNodeIdx();
    if(v->hasParent()){
      auto& parentsOutgroups = m_lcaOutgroups[v->getParent().lock()->getNodeIdx()];
      m_lcaOutgroups[vid].insert(parentsOutgroups.begin(),
                                 parentsOutgroups.end());
     // std::cout << "Np---" << v->getValue() << "p" << v->getParent().lock()->getValue();
     // for(auto species : m_lcaOutgroups[vid]){
     //   std::cout << "  " << species;
     // }
     // std::cout << std::endl;
    }

    if(v->getChildren().size() == 1){
      continue;
    }

    for(auto& c1 : v->getChildren()){
      size_t c1id = c1->getNodeIdx();
      for(auto& c2 : v->getChildren()){
        if(c1 == c2){
          continue;
        }
        size_t c2id = c2->getNodeIdx();

        if(rootDone){
          auto& genesC1 = m_subtreeGenes[c1id];
          auto& speciesC2 = m_subtreeSpecies[c2id];
          auto outgroupCandidates = std::unordered_set<std::shared_ptr<TreeNode>>(speciesC2.begin(),
                                                                                  speciesC2.end());
          if((genesC1.size() > 1) && (speciesC2.size() > 1)){
            for(size_t i = 0; i < speciesC2.size()-1; ++i){
              for(size_t j = i; j < speciesC2.size(); ++j){
                auto genesC2 = std::vector<Gene*>();
                auto& genesToAppend1 = m_ptrS->getSpeciesGenes(speciesC2[i]->getValue());
                genesC2.insert(genesC2.end(),
                               genesToAppend1.begin(),
                               genesToAppend1.end());
                auto& genesToAppend2 = m_ptrS->getSpeciesGenes(speciesC2[j]->getValue());
                genesC2.insert(genesC2.end(),
                               genesToAppend2.begin(),
                               genesToAppend2.end());

                std::vector<double> votes {0.0, 0.0};
//                for(size_t k = 0; k < m_ptrParam->getOutgroupLimit(); ++k){
                for(size_t k = 0; k < 20; ++k){

                  // sampling of a, b, c, d
                  if(m_benchmark) m_benchmark->startOutgroupsCorrectSample();
                  auto a_b = std::vector<Gene*>();
                  std::sample(genesC1.begin(), genesC1.end(),
                              std::back_inserter(a_b),
                              2, std::mt19937{std::random_device{}()});
                  auto c_d = std::vector<Gene*>();
                  std::sample(genesC2.begin(), genesC2.end(),
                              std::back_inserter(c_d),
                              2, std::mt19937{std::random_device{}()});
                  if(m_benchmark) m_benchmark->endOutgroupsCorrectSample();

                  if(m_ptrParam->weightedMode()){
                    auto voteAndWeight = m_ptrQ->supportedQuartetWeighted(a_b[0], a_b[1], c_d[0], c_d[1]);
                    if(voteAndWeight.first == 0){
                      votes[0] += voteAndWeight.second;
                    } else {
                      votes[1] += voteAndWeight.second;
                    }
                  } else {
                    if(m_ptrQ->supportedQuartetMajority(a_b[0], a_b[1], c_d[0], c_d[1]) == 0){
                      votes[0] += 1.0;
                    } else {
                      votes[1] += 1.0;
                    }
                  }
                }

                if(votes[1] / (votes[0]+votes[1]) >= m_ptrParam->getIncongruenceThreshold()){
                  outgroupCandidates.erase(speciesC2[i]);
                  outgroupCandidates.erase(speciesC2[j]);
                }
              }
            }
          }
          for(auto species : outgroupCandidates){
            m_lcaOutgroups[c1id].insert(species->getValue());
          }

        } else {
          // std::cout << "N---" << v->getValue();
          // for(auto species : m_subtreeSpecies[c2id]){
          //   std::cout << "  " << species->getValue();
          // }
          // std::cout << std::endl;
          for(auto species : m_subtreeSpecies[c2id]){
            m_lcaOutgroups[c1id].insert(species->getValue());
          }
        }

      }
    }
    rootDone = true;
  }
}

void
OutgroupChoice::initialize(){
  if(m_benchmark) m_benchmark->startSortDistances();
  buildIMatrix();
  if(m_benchmark) m_benchmark->endSortDistances();

  if(m_benchmark) {
    m_benchmark->startOutgroupInit();
    m_benchmark->startOutgroupsLca();
  }

  // compute the last common ancestors in S
  computeLcaS();

  if(m_benchmark) {
    m_benchmark->endOutgroupsLca();
    m_benchmark->startOutgroupsCorrect();
  }

  // assign corrected outgroups
  computeOutgroups();

  if(m_benchmark) {
    m_benchmark->endOutgroupsCorrect();
    m_benchmark->endOutgroupInit();
  }
}

std::vector<Gene*>
OutgroupChoice::getClosest(Gene* x, std::vector<Gene*>& genesY){

  size_t leafIdx1 = m_speciesLeafIdx[x->getSpecies()];
  size_t leafIdx2 = m_speciesLeafIdx[genesY[0]->getSpecies()];
  std::unordered_set<std::string>& outgroupsS = m_lcaOutgroups[m_lcaS.at(leafIdx1,leafIdx2)];

  auto result = std::vector<Gene*>();
  size_t i = x->getIndex();
  size_t j = 0;
  size_t N = m_I.getDim();
  size_t limit = m_ptrParam->getOutgroupLimit();

  while(result.size() < limit  && j < N){
  //while(j < N){
    Gene* candidateGene = m_ptrS->getGenePtr(m_I.at(i,j));
    if(outgroupsS.find(candidateGene->getSpecies()) != outgroupsS.end()) {
      result.push_back(candidateGene);
    }
    ++j;
  }
  // std::cout << x->getIdentifier()
  //           << "  " << m_lcaOutgroups[m_lcaS.at(leafIdx1,leafIdx2)].size()
  //           << "  " << genesY[0]->getSpecies()
  //           << "  " << m_nodeIdxToNode[m_lcaS.at(leafIdx1,leafIdx2)]
  //           << std::endl;
  // std::cout << result.size() << std::endl;
  // auto result2 = std::vector<Gene*>();
  // std::sample(result.begin(), result.end(),
  //             std::back_inserter(result2),
  //             std::min(m_ptrParam->getOutgroupLimit(), result.size()),
  //             std::mt19937{std::random_device{}()});
  //
  // return result2;
  return result;
}
