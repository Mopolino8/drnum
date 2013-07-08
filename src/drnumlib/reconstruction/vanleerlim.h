#ifndef VANLEERLIM_H
#define VANLEERLIM_H

#include "blockcfd.h"

struct VanLeerLim
{
  static CUDA_DH real lim(real delta1, real delta2)
  {
    countFlops(5);
    real r   = delta2/nonZero(delta1, GLOBAL_EPS);
    return (r + fabs(r))/(1 + fabs(r));
  }
};

#endif // VANLEERLIM_H
