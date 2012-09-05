#ifndef VANALBADA_H
#define VANALBADA_H

#include "blockcfd.h"

struct VanAlbada
{
  static real lim(real delta1, real delta2)
  {
    real r   = delta2/nonZero(delta1, global_eps);
    return (sqr(r) + r)/(sqr(r) + 1);
    //return max(real(0), (2*delta1*delta2/max(real(1e-3), sqr(delta1)+sqr(delta2))));
  }
};

#endif // VANALBADA_H
