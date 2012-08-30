#ifndef MINMOD_H
#define MINMOD_H

#include "blockcfd.h"

struct MinMod
{
  static real delta(real delta1, real delta2)
  {
    countFlops(1);
    if (delta1*delta2 <= 0) {
      return 0;
    } else if (fabs(delta1) > fabs(delta2)) {
      return delta2;
    }
    return delta1;
  }
};

#endif // MINMOD_H
