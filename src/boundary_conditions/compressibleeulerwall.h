#ifndef COMPRESSIBLEEULERWALL_H
#define COMPRESSIBLEEULERWALL_H

struct CompressibleEulerWall
{
  static void correct(real nx, real ny, real nz, real* var);
};

void CompressibleEulerWall::correct(real nx, real ny, real nz, real *var)
{
  var[1] -= var[1]*nx;
  var[2] -= var[2]*ny;
  var[3] -= var[3]*nz;
}

#endif // COMPRESSIBLEEULERWALL_H
