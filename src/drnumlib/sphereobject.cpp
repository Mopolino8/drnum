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
#include "sphereobject.h"


SphereObject::SphereObject() :
  ObjectDefinition()
{
}


void SphereObject::setParams (real xo_center, real yo_center, real zo_center, real radius)
{
  m_XoCenter = xo_center;
  m_YoCenter = yo_center;
  m_ZoCenter = zo_center;
  m_Radius = radius;
  m_QRadius = radius * radius;
}


bool SphereObject::isInside (const real &xo, const real &yo, const real &zo)
{
  // all on xyzo-coords

  real dx = xo - m_XoCenter;
  real dy = yo - m_YoCenter;
  real dz = zo - m_ZoCenter;
  real q_dist = dx*dx + dy*dy + dz*dz;

  return (q_dist < m_QRadius);
}
