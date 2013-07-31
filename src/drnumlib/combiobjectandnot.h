#ifndef COMBIOBJECTANDNOT_H
#define COMBIOBJECTANDNOT_H

class CombiObjectAndNot;

#include "combiobject.h"

class CombiObjectAndNot : public CombiObject
{
public:
  CombiObjectAndNot(ObjectDefinition* object_a, ObjectDefinition* object_b);
  CombiObjectAndNot(ObjectDefinition* object_a);
  virtual bool evalBool();
};

#endif // COMBIOBJECTANDNOT_H
