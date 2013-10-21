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
#include "combiobjectandnot.h"

CombiObjectAndNot::CombiObjectAndNot(ObjectDefinition* object_a, ObjectDefinition* object_b)
  : CombiObject(object_a, object_b)
{
}


CombiObjectAndNot::CombiObjectAndNot(ObjectDefinition* object_a)
  : CombiObject(object_a)
{
}


bool CombiObjectAndNot::evalBool()
{
  // Ask for: "the first object inside and not any other".

  bool inside = m_Objects[0]->evalBool();

  if (inside) {
    for (size_t i_o = 1; i_o < m_Objects.size(); i_o++) {
      if(m_Objects[i_o]->evalBool()) {
        inside = false;
        break;
      }
    }
  }

  return inside;


//  bool a_inside = m_ObjectA->evalBool();
//  bool b_inside = m_ObjectB->evalBool();

//  return a_inside && !b_inside;
}
