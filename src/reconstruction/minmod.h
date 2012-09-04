#ifndef MINMOD_H
#define MINMOD_H

#include "blockcfd.h"

struct MinMod
{
  static real lim(real delta1, real delta2)
  {
    return max(real(0), min(real(1), real(1)*delta2/nonZero(delta1, global_eps)));
  }
};

#endif // MINMOD_H
