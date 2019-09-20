#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

class Parameters {
public:
  Parameters()
   : m_matrixFile("")
   , m_speciesFile("")
   , m_treeFile("")
   , m_benchmarkFile("benchmark.txt")
   , m_disableQuartet(false)
   , m_restrictedY(false)
   , m_weightedMode(false)
   , m_subtreeFiles(false)
   , m_relativeOutgroups(false)
   , m_benchmark(false)
   , m_outgroupLimit(10)
   , m_epsilon(0.5)
   , m_incongruenceThreshold(0.2) { };

   void
   parseParameters(int argc, char* argv[]);

   bool
   checkIntegrity();

   const std::string&
   getMatrixFilename() const { return m_matrixFile; };

   const std::string&
   getSpeciesFilename() const { return m_speciesFile; };

   const std::string&
   getTreeFilename() const { return m_treeFile; };

   const std::string&
   getBenchmarkFilename() const { return m_benchmarkFile; };

   bool
   quartetsDisabled() const { return m_disableQuartet; };

   bool
   restrictedY() const { return m_restrictedY; };

   bool
   weightedMode() const { return m_weightedMode; };

   bool
   subtreeFiles() const { return m_subtreeFiles; };

   bool
   relativeOutgroups() const { return m_relativeOutgroups; };

   bool
   benchmark() const { return m_benchmark; };

   size_t
   getOutgroupLimit() const { return m_outgroupLimit; };

   double
   getEpsilon() const { return m_epsilon; };

   double
   getIncongruenceThreshold() const { return m_incongruenceThreshold; };

private:
  std::string m_matrixFile, m_speciesFile, m_treeFile, m_benchmarkFile;

  bool m_disableQuartet, m_restrictedY, m_weightedMode, m_subtreeFiles, m_relativeOutgroups, m_benchmark;
  size_t m_outgroupLimit;
  double m_epsilon, m_incongruenceThreshold;
};

#endif /* PARAMETERS_H */
