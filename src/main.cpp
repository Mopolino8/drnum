#ifndef GPU
//#include "examples/cpujet.h"
//#include "examples/cpujet_mb_grid1.h"
#include "examples/cpujet_mb_onlygrid.h"
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

