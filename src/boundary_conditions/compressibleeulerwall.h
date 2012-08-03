#ifndef COMPRESSIBLEEULERWALL_H
#define COMPRESSIBLEEULERWALL_H

#include "blockcfd.h"

struct CompressibleEulerWall
{
  static void correct(real nx, real ny, real nz, real* var);
};

void CompressibleEulerWall::correct(real nx, real ny, real nz, real *var)
{
  real Un = var[1]*nx + var[2]*ny + var[3]*nz;
  var[1] -= Un*nx;
  var[2] -= Un*ny;
  var[3] -= Un*nx;
  countFlops(11);
}

#endif // COMPRESSIBLEEULERWALL_H
