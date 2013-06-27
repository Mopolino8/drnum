#ifdef GPU
  #include "external_aero.h"
#endif

extern "C"
void GPU_main()
{
#ifdef GPU
  run();
#endif
}

