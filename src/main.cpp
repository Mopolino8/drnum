#ifndef GPU

// original
//#include "examples/cpujet.h"

// single patch on a PatchGrid
//#include "examples/cpujet_pg_single.h"

// two patches on a PatchGrid
//#include "examples/cpujet_pg_dual_1.h"

// everything from grid file, assuming iterators match
#include "examples/cpujet_pg.h"  // next to test

// ??
//#include "examples/cpujet_mb_grid1.h"

// PatchGrid with two blocks, no computing
//#include "examples/cpujet_mb_onlygrid.h"

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

