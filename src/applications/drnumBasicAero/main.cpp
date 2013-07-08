#ifndef GPU
#include "main.h"
#endif

extern "C" void GPU_main();

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

