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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef OBJECTDEFINITION_H
#define OBJECTDEFINITION_H

class ObjectDefinition;

#include "blockcfd.h"

class ObjectDefinition
{
  /** @todo There seams to be a bug in the C++ compiler. It will not
   *  accept other instances of derived classes of ObjectDefinition
   *  access the below functions, if in protected range (???) */

  //C++ bug? friend class ObjectDefinition;

protected: // attributes
  bool m_KnownInside;

protected: // methods
  //C++ bug? setKnownInside(const bool& inside) {m_KnownInside = inside;}

  /** @todo define wether to use o-coords, or generalized parent coords xyz
    *       in all derived class. */

public:

  void setKnownInside(const bool& inside) {m_KnownInside = inside;}
  bool getKnownInside() {return m_KnownInside;}

  ObjectDefinition();

  virtual bool isInside(const real& xo, const real& yo, const real& zo) {
    cout << " need derived class for ObjectDefinition::isInside" << endl;
    BUG;
  }

  bool isInside(vec3_t xyzo) { return isInside (xyzo[0], xyzo[1],xyzo[2]); }

  virtual void getLowestObjects(vector<ObjectDefinition*>& my_lowest_objects);

  virtual bool evalBool() {return getKnownInside();}

};

#endif // OBJECTDEFINITION_H
