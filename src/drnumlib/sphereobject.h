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
#ifndef SPHEREOBJECT_H
#define SPHEREOBJECT_H

class SphereObject;

#include "blockobject.h"

class SphereObject : public BlockObject
{

protected:
  real m_XoCenter;
  real m_YoCenter;
  real m_ZoCenter;
  real m_Radius;
  real m_QRadius;

public:

  SphereObject(PatchGrid* patch_grid);

  void setParams (real xo_center, real yo_center, real zo_center, real radius);

  virtual bool isInside (const real& xo, const real& yo, const real& zo);

};

#endif // SPHEREOBJECT_H
