#ifndef QUARTETS_H
#define QUARTETS_H

#include <vector>
#include <utility>

#include "Gene.h"
#include "Scenario.h"
#include "Parameters.h"

class Quartets {
public:
  Quartets(Scenario* ptrS,
           Parameters* ptrParam)
   : m_ptrS(ptrS)
   , m_ptrParam(ptrParam) { };

   std::vector<Gene*>
   findBestMatches(const Gene* x, const std::vector<Gene*>& genesY, const std::vector<Gene*>& outgroupsZ);

   size_t
   supportedQuartetMajority(const Gene* x, const Gene* y1, const Gene* y2, const Gene* z);

   std::pair<size_t,double>
   supportedQuartetWeighted(const Gene* x, const Gene* y1, const Gene* y2, const Gene* z);

private:
  Scenario* m_ptrS;
  Parameters* m_ptrParam;
};

#endif /* QUARTETS_H */
