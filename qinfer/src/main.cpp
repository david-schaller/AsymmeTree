#include <fstream>

#include "Scenario.h"
#include "BMGBuilder.h"
#include "Benchmark.h"
#include "Parameters.h"

auto
main(int argc, char* argv[]) -> int
{
  Parameters param;
  param.parseParameters(argc, argv);
  if(!param.checkIntegrity()) return -1;

  Benchmark bm(&param);
  bool benchmark(param.benchmark());

  // read the input files and check integrity
  if(benchmark) bm.startReadFiles();
  auto s = Scenario(&param, benchmark ? &bm : nullptr);
  s.parseFiles();
  if(benchmark) bm.endReadFiles();

  // initialize intance of BMGBuilder
  auto bmgBuilder = BMGBuilder(&s, &param, benchmark ? &bm : nullptr);
  bmgBuilder.buildBMG();
  bmgBuilder.printBMG();

  // write benchmarking results into file
  if(benchmark){
    bm.endTotal();
    std::ofstream ostrm(param.getBenchmarkFilename(), std::ios::app);
    bm.flush(ostrm, s);
  }

  return 0;
}
