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
#ifndef CYLINDERINCARTESIANPATCH_H
#define CYLINDERINCARTESIANPATCH_H

#include "genericoperation.h"
#include "cartesianpatch.h"

class CylinderInCartesianPatch : public GenericOperation
{

  CartesianPatch *m_Patch;
  real            m_Temp;
  real            m_Rad;
  real            m_Omega;
  vec3_t          m_XOrg;

  real   norm(vec3_t x);
  real   dot(vec3_t x1 , vec3_t x2);
  real   scalarProduct(vec3_t x1 , vec3_t x2);
  vec3_t normalize(vec3_t x);

public:

  /**
   * This class sets two layers outside a cylinder of diameter m_Rad to a given pressure and velocity.
   * For storage during the computation of the layer data the residual field is abused.
   * The 5 storage indices per point are used as followed:
   *  - 1 cell type: 1-> in cylinder, 2-> layer 1, 3-> layer 2, 0-> outside, set to zero.
   *  - 2 total weight, w1 + w2 + w3.., accounting for the contribution of the cylinder edge cells to the "bc"
   *  - 3 weighted pressure, w1*p1 + w2*p2..
   *  - 4 u
   *  - 5 v
   */

  CylinderInCartesianPatch(CartesianPatch* patch);

  virtual void operator()();

};


inline real CylinderInCartesianPatch::norm(vec3_t x)
{
  return sqrt(x[0]*x[0] + x[1]*x[1] +x[2]*x[2]);
}

inline vec3_t CylinderInCartesianPatch::normalize(vec3_t x)
{
  x[0] /= norm(x);
  x[1] /= norm(x);
  x[2] /= norm(x);
  return  x;
}

inline real CylinderInCartesianPatch::dot(vec3_t x1, vec3_t x2)
{
  return (x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]);
}

inline real CylinderInCartesianPatch::scalarProduct(vec3_t x1, vec3_t x2) {
  return CylinderInCartesianPatch::dot(x1, x2)/CylinderInCartesianPatch::norm(x2);
}

#endif // CYLINDERINCARTESIANPATCH_H
