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

/**
 * Class for outgroup choice.
 *
 * This class contains functions for the selection of the closest outgroups corrected with evidence for ancient
 * duplication events.
 */
class OutgroupChoice {
public:
  /**
   * Class constructor which creates a new instance.
   * @param ptrS pointer to the Scenario instance.
   * @param ptrParam pointer to the Parameters instance.
   * @param bm pointer to the Benchmark instance.
   */
  OutgroupChoice(Scenario* ptrS, Quartets* ptrQ,
                 Parameters* ptrParam,
                 Benchmark* bm = nullptr)
    : m_ptrS(ptrS)
    , m_ptrQ(ptrQ)
    , m_ptrParam(ptrParam)
    , m_benchmark(bm) { };

  /**
   * Initialization steps for the choice of the closest outgroups including the computation of last common ancestors
   * in the species tree S, the correction for ancient duplications and the sorting of distances for each gene x.
   */
  void
  initialize();

  /**
   * Returns the closest outgroups for a gene x.
   */
  std::vector<Gene*>
  getClosest(Gene* x, std::vector<Gene*>& genesY);

private:
  Scenario* m_ptrS;           /*!< pointer to the Scenario instance */
  Quartets* m_ptrQ;           /*!< pointer to the Quartets instance */
  Parameters* m_ptrParam;     /*!< pointer to the Parameters instance */
  Benchmark* m_benchmark;     /*!< pointer to the Benchmark instance */

  Matrix<size_t> m_I;         /*!< index matrix for the genes sorted by distance */
  Matrix<size_t> m_lcaS;      /*!< matrix for the last common ancestors of the species */

  std::unordered_map<std::string, size_t> m_speciesLeafIdx;             /*!< indeces of the species */
  std::vector<std::vector<std::shared_ptr<TreeNode>>> m_subtreeSpecies; /*!< species nodes in the subtrees of the nodes */
  std::vector<std::vector<Gene*>> m_subtreeGenes;                       /*!< genes in the (species) subtrees */
  std::unordered_map<size_t, std::unordered_set<std::string>> m_lcaOutgroups; /*!< valid outgroups per lca  */

  /**
   * Sorts a row of the index matrix based on the gene distances.
   */
  void quicksortRow(size_t row, int l, int r);

  /**
   * Initializes and sorts the index matrix (by the gene distances).
   */
  void buildIMatrix();

  /**
   * Computes the last common ancestors of the leaves of the species tree.
   */
  void computeLcaS();

  /**
   * Computes the possible outgroups for the internal nodes of the species tree with correction for ancient
   * duplication events.
   */
  void computeOutgroups();
};


#endif /* OUTGROUPCHOICE_H */
