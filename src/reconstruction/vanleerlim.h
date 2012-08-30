#ifndef VANLEERLIM_H
#define VANLEERLIM_H

#include "blockcfd.h"

struct VanLeerLim
{
  static real delta(real delta1, real delta2)
  {
    //return 2/(1/nonZero(delta1, global_eps) + 1/nonZero(delta2, global_eps));
    if (delta1*delta2 <= 0) {
      return 0;
    }
    countFlops(4);
    return 2/(1/delta1 + 1/delta2);
  }
};

#endif // VANLEERLIM_H
