#include "coordtransform.h"

CoordTransform::CoordTransform()
{
  axx = 1; axy = 0; axz = 0, bx = 0;
  ayx = 0; ayy = 1; ayz = 0, by = 0;
  azx = 0; azy = 0; azz = 1, bz = 0;
}

void CoordTransform::setMatrix(mat3_t A)
{
  axx = A[0][0]; axy = A[0][1]; axz = A[0][2];
  ayx = A[1][0]; ayy = A[1][1]; ayz = A[1][2];
  azx = A[2][0]; azy = A[2][1]; azz = A[2][2];
}

void CoordTransform::setVector(vec3_t b)
{
  bx = b[0]; by = b[1]; bz = b[2];
}

void CoordTransform::setAll(mat3_t A, vec3_t b)
{
  axx = A[0][0]; axy = A[0][1]; axz = A[0][2];
  ayx = A[1][0]; ayy = A[1][1]; ayz = A[1][2];
  azx = A[2][0]; azy = A[2][1]; azz = A[2][2];
  bx = b[0]; by = b[1]; bz = b[2];
}

void CoordTransform::setTransFromTo(const CoordTransform& c_from, const CoordTransform& c_to)
{
  copy((c_from.inverse()).concatenate(c_to));
}

mat3_t CoordTransform::getMatrix()
{ 
  mat3_t A;
  A[0][0] = axx; A[0][1] = axy; A[0][2] = axz;
  A[1][0] = ayx; A[1][1] = ayy; A[1][2] = ayz;
  A[2][0] = azx; A[2][1] = azy; A[2][2] = azz;
  return A;
}

vec3_t CoordTransform::getVector()
{
  return vec3_t(bx, by, bz);
}

CoordTransform CoordTransform::inverse() const
{
  CoordTransform t;
  mat3_t M;
  M[0][0] = axx; M[0][1] = axy; M[0][2] = axz;
  M[1][0] = ayx; M[1][1] = ayy; M[1][2] = ayz;
  M[2][0] = azx; M[2][1] = azy; M[2][2] = azz;
  M = M.inverse();
  t.axx = M[0][0]; t.axy = M[0][1]; t.axz = M[0][2];
  t.ayx = M[1][0]; t.ayy = M[1][1]; t.ayz = M[1][2];
  t.azx = M[2][0]; t.azy = M[2][1]; t.azz = M[2][2];
  vec3_t v(-bx, -by, -bz);
  vec3_t tv = M*v;
  t.bx = tv[0]; t.by = tv[1]; t.bz = tv[2];
  return t;
}
