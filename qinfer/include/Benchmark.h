#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <chrono>
#include <iosfwd>

#include "Parameters.h"

// forward declaration to avoid cyclic imports
class Scenario;

/**
 * Class for benchmarking.
 *
 * This class holds time points for benchmarking and contains functions for measuring the execution time.
 */
class Benchmark {
public:
  /**
   * Class constructor which creates a new instance.
   * @param ptrParam pointer to the Parameters instance.
   */
  Benchmark(Parameters* ptrParam)
   : m_ptrParam(ptrParam)
   , m_initTime(std::chrono::high_resolution_clock::now())
   , m_readFilesStart()
   , m_sortDistancesStart()
   , m_epsilonStart()
   , m_outgroupInitStart()
   , m_outgroupsLcaStart()
   , m_outgroupsCorrectStart()
   , m_outgroupsCorrectSampleStart()
   , m_buildBMGStart()
   , m_chooseOutrgroupsStart()

   , m_readFilesDuration()
   , m_sortDistancesDuration()
   , m_epsilonDuration()
   , m_outgroupInitDuration()
   , m_outgroupsLcaDuration()
   , m_outgroupsCorrectDuration()
   , m_outgroupsCorrectSampleDuration()
   , m_buildBMGDuration()
   , m_chooseOutrgroupsDuration() {
     // initialize start times to now, duration to 0 (this is a bit dirty)
     m_readFilesStart = std::chrono::high_resolution_clock::now();
     m_sortDistancesStart = std::chrono::high_resolution_clock::now();
     m_epsilonStart = std::chrono::high_resolution_clock::now();
     m_outgroupInitStart = std::chrono::high_resolution_clock::now();
     m_outgroupsLcaStart = std::chrono::high_resolution_clock::now();
     m_outgroupsCorrectStart = std::chrono::high_resolution_clock::now();
     m_outgroupsCorrectSampleStart = std::chrono::high_resolution_clock::now();
     m_buildBMGStart = std::chrono::high_resolution_clock::now();
     m_chooseOutrgroupsStart = std::chrono::high_resolution_clock::now();

     m_readFilesDuration = m_readFilesStart - m_readFilesStart;
     m_sortDistancesDuration = m_readFilesStart - m_readFilesStart;
     m_epsilonDuration = m_readFilesStart - m_readFilesStart;
     m_outgroupInitDuration = m_readFilesStart - m_readFilesStart;
     m_outgroupsLcaDuration = m_readFilesStart - m_readFilesStart;
     m_outgroupsCorrectDuration = m_readFilesStart - m_readFilesStart;
     m_outgroupsCorrectSampleDuration = m_readFilesStart - m_readFilesStart;
     m_buildBMGDuration = m_readFilesStart - m_readFilesStart;
     m_chooseOutrgroupsDuration = m_readFilesStart - m_readFilesStart;
     m_totalDuration = m_readFilesStart - m_readFilesStart;
    };

  /**
   * Start and end the time measurement for reading the input files.
   */
  void startReadFiles() { m_readFilesStart = std::chrono::high_resolution_clock::now(); };
  void endReadFiles() { m_readFilesDuration = std::chrono::high_resolution_clock::now() - m_readFilesStart; };

  /**
   * Start and end the time measurement for sorting the distances per gene.
   */
  void startSortDistances() { m_sortDistancesStart = std::chrono::high_resolution_clock::now(); };
  void endSortDistances() { m_sortDistancesDuration = std::chrono::high_resolution_clock::now() - m_sortDistancesStart; };

  /**
   * Start and end the time measurement of the Extended Best Hits method.
   */
  void startEpsilon() { m_epsilonStart = std::chrono::high_resolution_clock::now(); };
  void endEpsilon() { m_epsilonDuration = std::chrono::high_resolution_clock::now() - m_epsilonStart; };

  /**
   * Start and end the time measurement for the initialization of outgroups.
   */
  void startOutgroupInit() { m_outgroupInitStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupInit() { m_outgroupInitDuration = std::chrono::high_resolution_clock::now() - m_outgroupInitStart; };

  /**
   * Start and end the time measurement for computing last common ancestors in the species tree.
   */
  void startOutgroupsLca() { m_outgroupsLcaStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupsLca() { m_outgroupsLcaDuration = std::chrono::high_resolution_clock::now() - m_outgroupsLcaStart; };

  /**
   * Start and end the time measurement for the outgroup correction based on evidence for ancient duplications.
   */
  void startOutgroupsCorrect() { m_outgroupsCorrectStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupsCorrect() { m_outgroupsCorrectDuration = std::chrono::high_resolution_clock::now() - m_outgroupsCorrectStart; };

  /**
   * Start and end the time measurement for std::sample in outgroup correction (the time is cumulated).
   */
  void startOutgroupsCorrectSample() { m_outgroupsCorrectSampleStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupsCorrectSample() { m_outgroupsCorrectSampleDuration += std::chrono::high_resolution_clock::now() - m_outgroupsCorrectSampleStart; };

  /**
   * Start and end the time measurement for the inference of the Best Match Graph.
   */
  void startBuildBMG() { m_buildBMGStart = std::chrono::high_resolution_clock::now(); };
  void endBuildBMG() { m_buildBMGDuration = std::chrono::high_resolution_clock::now() - m_buildBMGStart; };

  /**
   * Start and end the time measurement for choosing outgroups (time is cumulated).
   */
  void startChooseOutgroups() { m_chooseOutrgroupsStart = std::chrono::high_resolution_clock::now(); };
  void endChooseOutgroups() { m_chooseOutrgroupsDuration += std::chrono::high_resolution_clock::now() - m_chooseOutrgroupsStart; };

  /**
   * End the time measurement of the overall execution of the inference.
   */
  void endTotal() { m_totalDuration += std::chrono::high_resolution_clock::now() - m_initTime; };

  /**
   * Writes the benchmarking results into a files.
   */
  void flush(std::basic_ofstream<char>& stream, Scenario& ptrS);

private:
  Parameters* m_ptrParam;

  std::chrono::time_point<std::chrono::high_resolution_clock> m_initTime;

  std::chrono::time_point<std::chrono::high_resolution_clock> m_readFilesStart,
                                                              m_sortDistancesStart,
                                                              m_epsilonStart,
                                                              m_outgroupInitStart,
                                                              m_outgroupsLcaStart,
                                                              m_outgroupsCorrectStart,
                                                              m_outgroupsCorrectSampleStart,
                                                              m_buildBMGStart,
                                                              m_chooseOutrgroupsStart;

  std::chrono::duration<double> m_readFilesDuration,
                                m_sortDistancesDuration,
                                m_epsilonDuration,
                                m_outgroupInitDuration,
                                m_outgroupsLcaDuration,
                                m_outgroupsCorrectDuration,
                                m_outgroupsCorrectSampleDuration,
                                m_buildBMGDuration,
                                m_chooseOutrgroupsDuration,
                                m_totalDuration;
};

#endif /* BENCHMARK_H */
