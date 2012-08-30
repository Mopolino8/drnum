#ifndef VANALBADA_H
#define VANALBADA_H

#include "blockcfd.h"

struct VanAlbada
{
  static real delta(real delta1, real delta2)
  {
    return 2*sqr(delta1)*delta2/nonZero(sqr(delta1) + sqr(delta2), global_eps);
  }
};

#endif // VANALBADA_H
