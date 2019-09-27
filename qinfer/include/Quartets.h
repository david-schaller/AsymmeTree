#ifndef QUARTETS_H
#define QUARTETS_H

#include <vector>
#include <utility>

#include "Gene.h"
#include "Scenario.h"
#include "Parameters.h"

/**
 * Class for the inference of quartet relations.
 *
 * This class contains functions for the inference of quartet relations and best matches.
 */
class Quartets {
public:
  /**
   * Class constructor which creates a new instance.
   * @param ptrS pointer to the Scenario instance.
   * @param ptrParam pointer to the Parameters instance.
   */
  Quartets(Scenario* ptrS,
           Parameters* ptrParam)
   : m_ptrS(ptrS)
   , m_ptrParam(ptrParam) { };

   /**
    * Finds the best matches for gene x in the candidate set Y using the set of outgroup genes Z.
    * @param x reference gene x.
    * @param genesY set of best match candidates Y.
    * @param outgroupsZ set of outgroup genes Z.
    */
   std::vector<Gene*>
   findBestMatches(const Gene* x, const std::vector<Gene*>& genesY, const std::vector<Gene*>& outgroupsZ);

   /**
    * Finds the supported quartet for four genes using majority voting.
    * @param x gene x.
    * @param y1 gene y1.
    * @param y2 gene y2.
    * @param z gene z.
    */
   size_t
   supportedQuartetMajority(const Gene* x, const Gene* y1, const Gene* y2, const Gene* z);

   /**
    * Finds the supported quartet for four genes using weighted voting.
    * @param x gene x.
    * @param y1 gene y1.
    * @param y2 gene y2.
    * @param z gene z.
    */
   std::pair<size_t,double>
   supportedQuartetWeighted(const Gene* x, const Gene* y1, const Gene* y2, const Gene* z);

private:
  Scenario* m_ptrS;           /*!< pointer to the Scenario instance */
  Parameters* m_ptrParam;     /*!< pointer to the Parameters instance */
};

#endif /* QUARTETS_H */
