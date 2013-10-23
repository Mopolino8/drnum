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
#ifndef LSOBCOPERATOR_H
#define LSOBCOPERATOR_H

class LSOBCOperator;

#include "lslayerdataextrapol.h"

/** @todo Maybe this class is suitable to hold access(const LSLayerdataExtrapol&)
  *       later. Then introduce template DIM .  */

/** Operator base class for LevelSetobject boundary conditions. */
class LSOBCOperator
{
public:
  LSOBCOperator();
};

LSOBCOperator::LSOBCOperator() {}

#endif // LSOBCOPERATOR_H
