#ifndef VANALBADA_H
#define VANALBADA_H

#include "blockcfd.h"

struct VanAlbada
{
  static real lim(real delta1, real delta2);
};

inline real VanAlbada::lim(real delta1, real delta2)
{
  countFlops(4);
  real r = delta2/nonZero(delta1, global_eps);
  return (sqr(r) + r)/(sqr(r) + 1);
}

#endif // VANALBADA_H
