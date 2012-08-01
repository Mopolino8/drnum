#ifndef COMPRESSIBLEEULERWALL_H
#define COMPRESSIBLEEULERWALL_H

struct CompressibleEulerWall
{
  static void correct(real nx, real ny, real nz, real* var);
};

void CompressibleEulerWall::correct(real nx, real ny, real nz, real *var)
{
  real rhoUn = var[1]*nx + var[2]*ny + var[3]*nz;
  var[1] -= rhoUn*nx;
  var[2] -= rhoUn*ny;
  var[3] -= rhoUn*nz;
}

#endif // COMPRESSIBLEEULERWALL_H
