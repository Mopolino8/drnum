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
