#ifndef OUTGROUPCHOICE_H
#define OUTGROUPCHOICE_H

#include <memory>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>

#include "Gene.h"
#include "Scenario.h"
#include "Matrix.h"
#include "Tree.h"
#include "Quartets.h"
#include "Parameters.h"
#include "Benchmark.h"

class OutgroupChoice {
public:
  OutgroupChoice(Scenario* ptrS, Quartets* ptrQ,
                 Parameters* ptrParam,
                 Benchmark* bm = nullptr)
    : m_ptrS(ptrS)
    , m_ptrQ(ptrQ)
    , m_ptrParam(ptrParam)
    , m_benchmark(bm) { };

  void
  initialize();

  std::vector<Gene*>
  getClosest(Gene* x, std::vector<Gene*>& genesY);

private:
  Scenario* m_ptrS;
  Quartets* m_ptrQ;
  Parameters* m_ptrParam;
  Benchmark* m_benchmark;

  Matrix<size_t> m_I;
  Matrix<size_t> m_lcaS;

  std::unordered_map<std::string, size_t> m_speciesLeafIdx;
  // std::unordered_map<size_t, std::string> m_nodeIdxToNode;
  std::vector<std::vector<std::shared_ptr<TreeNode>>> m_subtreeSpecies;
  std::vector<std::vector<Gene*>> m_subtreeGenes;
  std::unordered_map<size_t, std::unordered_set<std::string>> m_lcaOutgroups;

  void quicksortRow(size_t row, int l, int r);
  void buildIMatrix();
  void computeLcaS();
  void computeOutgroups();
};


#endif /* OUTGROUPCHOICE_H */
