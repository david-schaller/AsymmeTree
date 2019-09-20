#ifndef BENCHMARK_H
#define BENCHMARK_H

#include <chrono>
#include <iosfwd>

#include "Parameters.h"

class Scenario;

class Benchmark {
public:
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

  void startReadFiles() { m_readFilesStart = std::chrono::high_resolution_clock::now(); };
  void endReadFiles() { m_readFilesDuration = std::chrono::high_resolution_clock::now() - m_readFilesStart; };

  void startSortDistances() { m_sortDistancesStart = std::chrono::high_resolution_clock::now(); };
  void endSortDistances() { m_sortDistancesDuration = std::chrono::high_resolution_clock::now() - m_sortDistancesStart; };

  void startEpsilon() { m_epsilonStart = std::chrono::high_resolution_clock::now(); };
  void endEpsilon() { m_epsilonDuration = std::chrono::high_resolution_clock::now() - m_epsilonStart; };

  void startOutgroupInit() { m_outgroupInitStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupInit() { m_outgroupInitDuration = std::chrono::high_resolution_clock::now() - m_outgroupInitStart; };

  void startOutgroupsLca() { m_outgroupsLcaStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupsLca() { m_outgroupsLcaDuration = std::chrono::high_resolution_clock::now() - m_outgroupsLcaStart; };

  void startOutgroupsCorrect() { m_outgroupsCorrectStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupsCorrect() { m_outgroupsCorrectDuration = std::chrono::high_resolution_clock::now() - m_outgroupsCorrectStart; };

  // functions for std::sample in outgroup correction, the time is cumulated
  void startOutgroupsCorrectSample() { m_outgroupsCorrectSampleStart = std::chrono::high_resolution_clock::now(); };
  void endOutgroupsCorrectSample() { m_outgroupsCorrectSampleDuration += std::chrono::high_resolution_clock::now() - m_outgroupsCorrectSampleStart; };

  void startBuildBMG() { m_buildBMGStart = std::chrono::high_resolution_clock::now(); };
  void endBuildBMG() { m_buildBMGDuration = std::chrono::high_resolution_clock::now() - m_buildBMGStart; };

  // functions for choosing outgroups are called multiple times, the time is cumulated
  void startChooseOutgroups() { m_chooseOutrgroupsStart = std::chrono::high_resolution_clock::now(); };
  void endChooseOutgroups() { m_chooseOutrgroupsDuration += std::chrono::high_resolution_clock::now() - m_chooseOutrgroupsStart; };

  void endTotal() { m_totalDuration += std::chrono::high_resolution_clock::now() - m_initTime; };

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
