#ifndef GPU

// original
//#include "examples/cpujet.h"

// single patch on a PatchGrid
#include "examples/cpujet_pg_single.h"  // next to test

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
  omp_set_num_threads(4);
  run();
#endif
}

