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
#ifndef COMBIOBJECT_H
#define COMBIOBJECT_H

class CombiObject;

#include "objectdefinition.h"

class CombiObject : public ObjectDefinition
{

protected: // attributes
  vector<ObjectDefinition*> m_LowestObjects;
  ObjectDefinition* m_ObjectA;
  ObjectDefinition* m_ObjectB;

  vector<ObjectDefinition*> m_Objects;


protected: // methods
  void considerLowestObjectsOf(ObjectDefinition* object);
  void concatLowestObjects(vector<ObjectDefinition*>& other_lowest_objects);
  void findLowestObjects();

public:
  CombiObject(ObjectDefinition* object_a, ObjectDefinition* object_b);
  CombiObject(ObjectDefinition* object_a);
  void includeObject(ObjectDefinition* object);
  virtual bool isInside(const real& xo, const real& yo, const real& zo);
  virtual bool evalBool() = 0;

};

#endif // COMBIOBJECT_H
