#ifdef GPU
#include "external_aero.h"
//#include "examples/gpujet.h"
#endif

extern "C"
void GPU_main()
{
#ifdef GPU
  run();
#endif
}

