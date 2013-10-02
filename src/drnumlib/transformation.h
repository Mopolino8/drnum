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
#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include "drnum.h"

class Transformation
{

  real axx, axy, axz, ax4;
  real ayx, ayy, ayz, ay4;
  real azx, azy, azz, az4;
  real a4x, a4y, a4z, a44;

public:

  Transformation();

  void setMatrix(mat3_t A);
  void setVector(vec3_t b);
  mat3_t getMatrix();
  vec3_t getVector();

  void operator()(real x, real y, real z, real& xt, real& yt, real& zt) const;
  void operator()(real& x, real& y, real& z) const { operator()(x, y, z, x, y, z); }
  vec3_t operator()(const vec3_t& x) const;
  Transformation combine(const Transformation& c) const;
  Transformation inverse() const;

};


inline void Transformation::operator()(real x, real y, real z, real& xt, real& yt, real& zt) const
{
  xt = axx*x + axy*y + axz*z + ax4;
  yt = ayx*x + ayy*y + ayz*z + ay4;
  zt = azx*x + azy*y + azz*z + az4;
  countFlops(18);
}

inline vec3_t Transformation::operator()(const vec3_t& x) const
{
  vec3_t xt;
  operator()(x[0], x[1], x[2], xt[0], xt[1], xt[2]);
  return xt;
}

inline Transformation Transformation::combine(const Transformation &c) const
{
  Transformation t;

  t.axx = axx*c.axx + axy*c.ayx + axz*c.azx + ax4*c.a4x;
  t.axy = axx*c.axy + axy*c.ayy + axz*c.azy + ax4*c.a4y;
  t.axz = axx*c.axz + axy*c.ayz + axz*c.azz + ax4*c.a4z;
  t.ax4 = axx*c.ax4 + axy*c.ay4 + axz*c.az4 + ax4*c.a44;
  countFlops(28);

  t.ayx = ayx*c.axx + ayy*c.ayx + ayz*c.azx + ay4*c.a4x;
  t.ayy = ayx*c.axy + ayy*c.ayy + ayz*c.azy + ay4*c.a4y;
  t.ayz = ayx*c.axz + ayy*c.ayz + ayz*c.azz + ay4*c.a4z;
  t.ay4 = ayx*c.ax4 + ayy*c.ay4 + ayz*c.az4 + ay4*c.a44;
  countFlops(28);

  t.azx = azx*c.axx + azy*c.ayx + azz*c.azx + az4*c.a4x;
  t.azy = azx*c.axy + azy*c.ayy + azz*c.azy + az4*c.a4y;
  t.azz = azx*c.axz + azy*c.ayz + azz*c.azz + az4*c.a4z;
  t.az4 = azx*c.ax4 + azy*c.ay4 + azz*c.az4 + az4*c.a44;
  countFlops(28);

  t.a4x = a4x*c.axx + a4y*c.ayx + a4z*c.azx + a44*c.a4x;
  t.a4y = a4x*c.axy + a4y*c.ayy + a4z*c.azy + a44*c.a4y;
  t.a4z = a4x*c.axz + a4y*c.ayz + a4z*c.azz + a44*c.a4z;
  t.a44 = a4x*c.ax4 + a4y*c.ay4 + a4z*c.az4 + a44*c.a44;
  countFlops(28);

  return t;
}

#endif // TRANSFORMATION_H
