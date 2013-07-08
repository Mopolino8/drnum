#ifdef GPU
#include "main.h"
#endif

extern "C"
void GPU_main()
{
#ifdef GPU
  run();
#endif
}

