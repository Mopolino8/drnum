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
#include "coneobject.h"

ConeObject::ConeObject() :
  ObjectDefinition()
{
}


void ConeObject::setParams(real xo_bottom, real yo_bottom, real zo_bottom,
                           real axis_xo, real axis_yo, real axis_zo,
                           real radius_bottom, real radius_top)
{
  m_BottomO[0] = xo_bottom;
  m_BottomO[1] = yo_bottom;
  m_BottomO[2] = zo_bottom;

  m_AxisO[0] = axis_xo;
  m_AxisO[1] = axis_yo;
  m_AxisO[2] = axis_zo;

  m_RadiusBottom = radius_bottom;
  m_RadiusTop = radius_top;

  m_QLength = m_AxisO.abs2();
  m_InvQLength = real(1.)/m_QLength;
  m_Length = sqrt(m_QLength);
  m_InvLength = real(1.)/m_Length;
}


bool ConeObject::isInside (const real &xo, const real &yo, const real &zo)
{
  // all on xyzo-coords

  vec3_t point_xyzo;
  point_xyzo[0] = xo;
  point_xyzo[1] = yo;
  point_xyzo[2] = zo;

  // check between bottom and top
  vec3_t bot2point = point_xyzo - m_BottomO;
  real scal = m_AxisO * bot2point;
  real scal_n = scal * m_InvQLength;
  if(scal_n < 0.) {
    return false;
  }
  else if(scal_n > 1.) {
    return false;
  }

  // check in radial limits
  //.. square of distance of point to axis
  vec3_t ax_cross_b2p = m_AxisO.cross(bot2point);
  real qdist_qaxilen = ax_cross_b2p.abs2();
  real qdist = qdist_qaxilen * m_InvQLength;
  //.. square of hull distance to axis at rel. coord
  real r_lim = (1. - scal_n) * m_RadiusBottom + scal_n * m_RadiusTop;
  real qr_lim = r_lim * r_lim;
  //.. check
  if (qdist > qr_lim) {
    return false;
  }

  // survived until here? => inside
  return true;
}
