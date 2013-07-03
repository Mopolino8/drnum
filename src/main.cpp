#ifndef GPU
// #include "external_aero.h"
#include "ghostfluid_test.h"
#endif

extern "C" void GPU_main();

//#include "examples/cpujet_mb_grid1.h"
//#include "examples/cpujet_mb_grid1.h"

int main()
{
#ifdef GPU
  GPU_main();
#else
  int num_threads = 2;
  omp_set_num_threads(num_threads);
  cout << endl;
  cout << "*** NUMBER THREADS: " << num_threads << endl;
  cout << endl;
  run();
#endif
}

