#ifndef ROELIM_H
#define ROELIM_H

#include "blockcfd.h"

struct RoeLim
{
  static CUDA_DH real delta(real delta1, real delta2)
  {
    countFlops(1);
    if (delta1*delta2 <= 0) {
      return 0;
    } else {
      countFlops(7);
      static const real eps = 1e-6;
      real D1  = nonZero(delta1, eps);
      real D2  = nonZero(delta2, eps);
      real lim = min(real(1), min(fabs(2*D1/(D1+D2)), fabs(2*D2/(D1+D2))));
      return lim*delta1;
    }
  }
};

#endif // ROELIM_H
