#ifndef COMPRESSIBLEEULERWALL_H
#define COMPRESSIBLEEULERWALL_H

#include "blockcfd.h"

struct CompressibleEulerWall
{
  static void correct(real nx, real ny, real nz, real* var);
};

void CompressibleEulerWall::correct(real nx, real ny, real nz, real *var)
{
  real Mom_n = var[1]*nx + var[2]*ny + var[3]*nz;
  var[1] -= Mom_n*nx;
  var[2] -= Mom_n*ny;
  var[3] -= Mom_n*nz;
  //countFlops(11);
  var[4] -= 0.5*sqr(Mom_n)/var[0];
  countFlops(14);
}

#endif // COMPRESSIBLEEULERWALL_H
