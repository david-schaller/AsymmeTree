#include "Quartets.h"

#include "DiGraph.h"
#include <cmath>

std::vector<Gene*>
Quartets::findBestMatches(const Gene* x,
                          const std::vector<Gene*>& genesY,
                          const std::vector<Gene*>& outgroupsZ){
  auto gamma = DiGraph<Gene*>();

  for(size_t i = 0; i < genesY.size(); ++i){
    Gene* y1 = genesY[i];
    for(size_t j = i+1; j < genesY.size(); ++j){
      Gene* y2 = genesY[j];

      std::vector<double> votes = {0.0, 0.0, 0.0, 0.0};

      for(const Gene* z : outgroupsZ){
        if(m_ptrParam->weightedMode()){
          auto voteAndWeight = supportedQuartetWeighted(x, y1, y2, z);
          votes[voteAndWeight.first] += voteAndWeight.second;
        } else {
          votes[supportedQuartetMajority(x, y1, y2, z)] += 1.0;
        }
      }

      int quartet = std::distance(votes.begin(),
                                  std::max_element(votes.begin(), votes.end()));
      switch(quartet){
        case 0:
          // 0: xy1 | y2z
          gamma.addEdge(y2, y1);
          break;
        case 1:
          // 1: xy2 | y1z
          gamma.addEdge(y1, y2);
          break;
        default:
          // 2: xz | y1y2   or   3: star
          gamma.addEdge(y2, y1);
          gamma.addEdge(y1, y2);
      }
    }
  }

  return gamma.getSccWithoutOutedges();
}

size_t
Quartets::supportedQuartetMajority(const Gene* x, const Gene* y1,
                                   const Gene* y2, const Gene* z){
  size_t xIdx = x->getIndex();
  size_t y1Idx = y1->getIndex();
  size_t y2Idx = y2->getIndex();
  size_t zIdx = z->getIndex();
  size_t quartet;

  double xy1_y2z = m_ptrS->getDistance(xIdx, y1Idx) + m_ptrS->getDistance(y2Idx, zIdx);
  double xy2_y1z = m_ptrS->getDistance(xIdx, y2Idx) + m_ptrS->getDistance(y1Idx, zIdx);
  double xz_y1y2 = m_ptrS->getDistance(xIdx, zIdx) + m_ptrS->getDistance(y1Idx, y2Idx);

  // 0: xy1 | y2z
  if(xy1_y2z < xy2_y1z && xy1_y2z < xz_y1y2){
    quartet = 0;
  // 1: xy2 | y1z
  } else if(xy2_y1z < xy1_y2z && xy2_y1z < xz_y1y2){
    quartet = 1;
  // 2: xz | y1y2
  } else if(xz_y1y2 < xy1_y2z && xz_y1y2 < xy2_y1z){
    quartet = 2;
  // 3: star topology
  } else {
    quartet = 3;
  }

  return quartet;
}

std::pair<size_t,double>
Quartets::supportedQuartetWeighted(const Gene* x, const Gene* y1,
                                   const Gene* y2, const Gene* z){
  size_t xIdx = x->getIndex();
  size_t y1Idx = y1->getIndex();
  size_t y2Idx = y2->getIndex();
  size_t zIdx = z->getIndex();
  size_t quartet;
  double weight = 0.0;

  double xy1_y2z = m_ptrS->getDistance(xIdx, y1Idx) + m_ptrS->getDistance(y2Idx, zIdx);
  double xy2_y1z = m_ptrS->getDistance(xIdx, y2Idx) + m_ptrS->getDistance(y1Idx, zIdx);
  double xz_y1y2 = m_ptrS->getDistance(xIdx, zIdx) + m_ptrS->getDistance(y1Idx, y2Idx);

  // 0: xy1 | y2z
  if(xy1_y2z < xy2_y1z && xy1_y2z < xz_y1y2){
    quartet = 0;
  // 1: xy2 | y1z
  } else if(xy2_y1z < xy1_y2z && xy2_y1z < xz_y1y2){
    quartet = 1;
  // 2: xz | y1y2
  } else if(xz_y1y2 < xy1_y2z && xz_y1y2 < xy2_y1z){
    quartet = 2;
  // 3: star topology
  } else {
    quartet = 3;
  }

  std::vector<double> sums = {xy1_y2z, xy2_y1z, xz_y1y2};
  std::sort(sums.begin(), sums.end());

  if(sums[2] > 0){
    weight = (1 - sums[0]/sums[2]) * std::exp(sums[1]-sums[2]);
  }

  return std::make_pair(quartet, weight);
}
