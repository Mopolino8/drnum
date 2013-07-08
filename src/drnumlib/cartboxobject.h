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
#ifndef CARTBOXOBJECT_H
#define CARTBOXOBJECT_H

class CartboxObject;

#include "blockobject.h"

class CartboxObject : public BlockObject
{

protected:
  real m_Xo_min;
  real m_Xo_max;
  real m_Yo_min;
  real m_Yo_max;
  real m_Zo_min;
  real m_Zo_max;

public:

  CartboxObject(PatchGrid* patch_grid);

  void setParams (real xo_min, real xo_max,
                  real yo_min, real yo_max,
                  real zo_min, real zo_max);

  virtual bool isInside (const real& xo, const real& yo, const real& zo);


};

#endif // CARTBOXOBJECT_H
