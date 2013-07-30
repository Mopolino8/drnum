#ifndef COMBIOBJECTANDNOT_H
#define COMBIOBJECTANDNOT_H

class CombiObjectAndNot;

#include "combiobject.h"

class CombiObjectAndNot : public CombiObject
{
public:
  CombiObjectAndNot(ObjectDefinition* object_a, ObjectDefinition* object_b);
  virtual bool evalBool();
};

#endif // COMBIOBJECTANDNOT_H
