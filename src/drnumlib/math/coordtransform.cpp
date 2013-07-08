// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "coordtransform.h"

CoordTransform::CoordTransform()
{
  axx = 1; axy = 0; axz = 0; bx = 0;
  ayx = 0; ayy = 1; ayz = 0; by = 0;
  azx = 0; azy = 0; azz = 1; bz = 0;
}

void CoordTransform::setMatrix(mat3_t A)
{
  axx = A[0][0]; axy = A[0][1]; axz = A[0][2];
  ayx = A[1][0]; ayy = A[1][1]; ayz = A[1][2];
  azx = A[2][0]; azy = A[2][1]; azz = A[2][2];
}

real CoordTransform::setMatrixFromBaseIJ(vec3_t base_i, vec3_t base_j)
{
  // normalize inputs base_i and base_j
  base_i.normalise();
  base_j.normalise();
  // compute base_k
  vec3_t base_k = base_i.cross(base_j);
  // length of base_k is measure for linear independence of base_i and base_j
  real base_k_len = base_k.abs();
  real health = base_k_len;
  // make base_j orthogonal to base_i
  real scal_ij = base_i * base_j;
  base_j = base_j - scal_ij * base_i;
  // re-normalize orthogonalised base_j and normalize base_k vectors (base_i unchanged)
  base_j.normalise();
  base_k = (1./base_k_len) * base_k;
  // fill in matrix
  axx = base_i[0]; axy = base_j[0]; axz = base_k[0];
  ayx = base_i[1]; ayy = base_j[1]; ayz = base_k[1];
  azx = base_i[2]; azy = base_j[2]; azz = base_k[2];
  return health;
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

void CoordTransform::scaleVector(real scfactor)
{
  bx *= scfactor;
  by *= scfactor;
  bz *= scfactor;
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
