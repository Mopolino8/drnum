#ifndef MINMOD_H
#define MINMOD_H

#include "blockcfd.h"

struct MinMod
{
  static CUDA_DH real lim(real delta1, real delta2)
  {
    return max(real(0), min(real(1), real(1)*delta2/nonZero(delta1, GLOBAL_EPS)));
  }
};

#endif // MINMOD_H
