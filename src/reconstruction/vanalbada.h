#ifndef VANALBADA_H
#define VANALBADA_H

#include "blockcfd.h"

struct VanAlbada
{
  static real delta(real delta1, real delta2)
  {
    return delta1*(2*delta1*delta2/max(real(1e-6), sqr(delta1)+sqr(delta2)));
  }
};

#endif // VANALBADA_H
