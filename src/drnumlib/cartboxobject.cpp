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
#include "cartboxobject.h"


CartboxObject::CartboxObject() :
  ObjectDefinition()
{
}


void CartboxObject::setParams (real xo_min, real xo_max,
                               real yo_min, real yo_max,
                               real zo_min, real zo_max)
{
  m_Xo_min = xo_min;
  m_Xo_max = xo_max;
  m_Yo_min = yo_min;
  m_Yo_max = yo_max;
  m_Zo_min = zo_min;
  m_Zo_max = zo_max;
}


bool CartboxObject::isInside (const real &xo, const real &yo, const real &zo)
{
  bool inside =
      (xo > m_Xo_min) &&
      (xo < m_Xo_max) &&
      (yo > m_Yo_min) &&
      (yo < m_Yo_max) &&
      (zo > m_Zo_min) &&
      (zo < m_Zo_max);

  return inside;
}
