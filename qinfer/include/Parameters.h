#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>

/**
 * Class for parameters.
 *
 * This class holds the user-specified parameters.
 */
class Parameters {
public:
  /**
   * Class constructor which creates a new instance. Initializes all parameters to the default values.
   */
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

   /**
    * Parses the parameters from the command line.
    * @param argc number of parameters.
    * @param argv array of parameters.
    */
   void
   parseParameters(int argc, char* argv[]);

   /**
    * Checks if parameters are not contradictory and the files' existence.
    */
   bool
   checkIntegrity();

   /**
    * Returns the name of the matrix file.
    */
   const std::string&
   getMatrixFilename() const { return m_matrixFile; };

   /**
    * Returns the name of the file containing the species identifiers with the corresponding genes.
    */
   const std::string&
   getSpeciesFilename() const { return m_speciesFile; };

   /**
    * Returns the name of the species tree file.
    */
   const std::string&
   getTreeFilename() const { return m_treeFile; };

   /**
    * Returns the name of the benchmark file.
    */
   const std::string&
   getBenchmarkFilename() const { return m_benchmarkFile; };

   /**
    * Returns true if the Quartet method is disabled (only Extended Best Hits method).
    */
   bool
   quartetsDisabled() const { return m_disableQuartet; };

   /**
    * Returns true if the best match candidates (Y) are restricted by a pre-filtering.
    */
   bool
   restrictedY() const { return m_restrictedY; };

   /**
    * Returns true if the weighted mode is enabled instead of the majority voting.
    */
   bool
   weightedMode() const { return m_weightedMode; };

   /**
    * Returns true if not a Newick tree is provided but a file containing the leaves in the subtrees of the root
    * of the species tree (per line).
    */
   bool
   subtreeFiles() const { return m_subtreeFiles; };

   /**
    * Returns true if all relative outgroups in the species tree are considered, and corrected for ancient
    * duplication events.
    */
   bool
   relativeOutgroups() const { return m_relativeOutgroups; };

   /**
    * Returns true if time benchmarking is enabled.
    */
   bool
   benchmark() const { return m_benchmark; };

   /**
    * Returns the limit for the number of outgroups (|Z|).
    */
   size_t
   getOutgroupLimit() const { return m_outgroupLimit; };

   /**
    * Returns the value of the parameter epsilon for the Extended Best Hits method.
    */
   double
   getEpsilon() const { return m_epsilon; };

   /**
    * Returns the value of the threshold parameter for discarding relative outgroup species.
    */
   double
   getIncongruenceThreshold() const { return m_incongruenceThreshold; };

private:
  std::string m_matrixFile, m_speciesFile, m_treeFile, m_benchmarkFile;

  bool m_disableQuartet, m_restrictedY, m_weightedMode, m_subtreeFiles, m_relativeOutgroups, m_benchmark;
  size_t m_outgroupLimit;
  double m_epsilon, m_incongruenceThreshold;
};

#endif /* PARAMETERS_H */
