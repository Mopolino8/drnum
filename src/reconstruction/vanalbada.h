#ifndef VANALBADA_H
#define VANALBADA_H

#include "blockcfd.h"

struct VanAlbada
{
  static CUDA_DH real lim(real delta1, real delta2)
  {
    countFlops(4);
    real r = delta2/nonZero(delta1, GLOBAL_EPS);
    return (sqr(r) + r)/(sqr(r) + 1);
  }
};

#endif // VANALBADA_H
