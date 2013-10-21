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
#include "cylinderobject.h"


CylinderObject::CylinderObject() :
  ObjectDefinition()
{
}


void CylinderObject::setParams(real xo_bottom, real yo_bottom, real zo_bottom,
                               real axis_xo, real axis_yo, real axis_zo,
                               real radius)
{
  m_BottomO[0] = xo_bottom;
  m_BottomO[1] = yo_bottom;
  m_BottomO[2] = zo_bottom;

  m_AxisO[0] = axis_xo;
  m_AxisO[1] = axis_yo;
  m_AxisO[2] = axis_zo;

  m_QLength = m_AxisO.abs2();

  m_Radius = radius;
  m_QRadius = radius * radius;

  m_QLR = m_QLength * m_QRadius;
}


bool CylinderObject::isInside (const real &xo, const real &yo, const real &zo)
{
  // all on xyzo-coords

  vec3_t point_xyzo;
  point_xyzo[0] = xo;
  point_xyzo[1] = yo;
  point_xyzo[2] = zo;

  // check between bottom and top
  vec3_t bot2point = point_xyzo - m_BottomO;
  real scal = m_AxisO * bot2point;
  if(scal < 0.) {
    return false;
  }
  else if(scal > m_QLength) {
    return false;
  }

  // check in radial limits
  vec3_t cross = m_AxisO.cross(bot2point);
  real q_dist_q_axilen = cross.abs2();  /// @todo this has unnecessary sqrt()
  if (q_dist_q_axilen > m_QLR) {
    return false;
  }

  // survived until here? => inside
  return true;
}
