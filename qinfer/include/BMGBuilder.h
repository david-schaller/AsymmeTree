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

/**
 * Class for BMG inference.
 *
 * This class contains functions for the inference of best matches from evolutionary distance data.
 */
class BMGBuilder {
public:
  /**
   * Class constructor which creates a new instance.
   * @param ptrS pointer to the Scenario instance.
   * @param ptrParam pointer to the Parameters instance.
   * @param bm pointer to the Benchmark instance.
   */
  BMGBuilder(Scenario* ptrS,
             Parameters* ptrParam,
             Benchmark* bm = nullptr)
   : m_ptrS(ptrS)
   , m_ptrParam(ptrParam)
   , m_quartets(Quartets(ptrS, ptrParam))
   , m_outgroupChoice (OutgroupChoice(ptrS, &m_quartets, ptrParam, bm))
   , m_benchmark(bm) { };

 /**
  * Builds a BMG from the distance data.
  */
  void buildBMG();

  /**
   * Prints the edges of the inferred BMG.
   */
  void printBMG();

private:
  Scenario* m_ptrS;                   /*!< pointer to the Scenario instance */
  Parameters* m_ptrParam;             /*!< pointer to the Parameters instance */
  Quartets m_quartets;                /*!< Quartets instance */
  OutgroupChoice m_outgroupChoice;    /*!< OutgroupChoice instance */
  Benchmark* m_benchmark;             /*!< pointer to the Benchmark instance */

  Matrix<int> m_bmCandidates;         /*!< matrix of best match candidates */
  DiGraph<Gene*> m_bmg;               /*!< Best Match Graph */

  /*********************************************************************************************************************
  *                                                COMMON FUCTIONS
  *                                         incl. Extended Best Hits method
  *********************************************************************************************************************/
  /**
   * Infers best matches using the Extended Best Hits method.
   */
  void
  epsilonMethod();

  /**
   * Pre-filtering step for the Quartets approach based on the Extended Best Hits method.
   */
  void
  buildCandidateMatrix();

  /*********************************************************************************************************************
                                                  OUTGROUG METHOD I
                                        (choose outgroups w.r.t the root of S)
  *********************************************************************************************************************/
  /**
   * Samples outgroups w.r.t the root of S.
   */
  std::vector<Gene*>
  chooseOutgroups(const std::vector<Gene*>& outgroupCandidates);

  /**
   * Infers best matches using the Quartet approach using randomly sampled outgroups w.r.t the root of S.
   */
  void
  buildWithRootOutgroups();

  /*********************************************************************************************************************
                                                 OUTGROUG METHOD II
                             (relative outgroups corrected with incongruent quartets)
  *********************************************************************************************************************/
  /**
   * Infers best matches using the Quartet approach using the closest (corrected outgroup) genes for each gene x.
   */
  void
  buildWithRelativeOutgroups();
};

#endif /* BMGBUILDER_H */
