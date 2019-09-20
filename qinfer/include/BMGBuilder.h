#ifndef BMGBUILDER_H
#define BMGBUILDER_H

#include <vector>
#include <utility>

#include "Gene.h"
#include "Scenario.h"
#include "DiGraph.h"
#include "Matrix.h"
#include "Quartets.h"
#include "OutgroupChoice.h"
#include "Parameters.h"
#include "Benchmark.h"

class BMGBuilder {
public:
  BMGBuilder(Scenario* ptrS,
             Parameters* ptrParam,
             Benchmark* bm = nullptr)
   : m_ptrS(ptrS)
   , m_ptrParam(ptrParam)
   , m_quartets(Quartets(ptrS, ptrParam))
   , m_outgroupChoice (OutgroupChoice(ptrS, &m_quartets, ptrParam, bm))
   , m_benchmark(bm) { };

  void buildBMG();
  void printBMG();

private:
  Scenario* m_ptrS;
  Parameters* m_ptrParam;
  Quartets m_quartets;
  OutgroupChoice m_outgroupChoice;
  Benchmark* m_benchmark;

  Matrix<int> m_bmCandidates;
  DiGraph<Gene*> m_bmg;

  // COMMON FUNCTIONS
  void
  epsilonMethod();

  void
  buildCandidateMatrix();

  // OUTGROUP METHOD I
  std::vector<Gene*>
  chooseOutgroups(const std::vector<Gene*>& outgroupCandidates);

  void
  buildRootOutgroups();

  // OUTGROUP METHOD II
  void
  buildRelativeOutgroups();
};

#endif /* BMGBUILDER_H */
