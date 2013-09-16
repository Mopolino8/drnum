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
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "transformation.h"

Transformation::Transformation()
{
  axx = 1; axy = 0; axz = 0, ax4 = 0;
  ayx = 0; ayy = 1; ayz = 0, ay4 = 0;
  azx = 0; azy = 0; azz = 1, az4 = 0;
  a4x = 0; a4y = 0; a4z = 0, a44 = 1;
}

void Transformation::setMatrix(mat3_t A)
{
  axx = A[0][0]; axy = A[0][1]; axz = A[0][2];
  ayx = A[1][0]; ayy = A[1][1]; ayz = A[1][2];
  azx = A[2][0]; azy = A[2][1]; azz = A[2][2];
}

void Transformation::setVector(vec3_t b)
{
  ax4 = b[0]; ay4 = b[1]; az4 = b[2];
}

mat3_t Transformation::getMatrix()
{
  mat3_t A;
  A[0][0] = axx; A[0][1] = axy; A[0][2] = axz;
  A[1][0] = ayx; A[1][1] = ayy; A[1][2] = ayz;
  A[2][0] = azx; A[2][1] = azy; A[2][2] = azz;
  return A;
}

vec3_t Transformation::getVector()
{
  return vec3_t(ax4, ay4, az4);
}

Transformation Transformation::inverse() const
{
  Transformation t;
  mat4_t M;
  M[0][0] = axx; M[0][1] = axy; M[0][2] = axz; M[0][3] = ax4;
  M[1][0] = ayx; M[1][1] = ayy; M[1][2] = ayz; M[1][3] = ay4;
  M[2][0] = azx; M[2][1] = azy; M[2][2] = azz; M[2][3] = az4;
  M[3][0] = a4x; M[3][1] = a4y; M[3][2] = a4z; M[3][3] = a44;
  M = M.inverse();
  t.axx = M[0][0]; t.axy = M[0][1]; t.axz = M[0][2], t.ax4 = M[0][3];
  t.ayx = M[1][0]; t.ayy = M[1][1]; t.ayz = M[1][2], t.ay4 = M[1][3];
  t.azx = M[2][0]; t.azy = M[2][1]; t.azz = M[2][2], t.az4 = M[2][3];
  t.a4x = M[3][0]; t.a4y = M[3][1]; t.a4z = M[3][2], t.a44 = M[3][3];
  return t;
}

