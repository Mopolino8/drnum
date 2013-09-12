#ifndef COMBIOBJECTOR_H
#define COMBIOBJECTOR_H

class CombiObjectOr;

#include "combiobject.h"

class CombiObjectOr : public CombiObject
{
public:
    CombiObjectOr(ObjectDefinition* object_a, ObjectDefinition* object_b);
    CombiObjectOr(ObjectDefinition* object_a);
    virtual bool evalBool();
};

#endif // COMBIOBJECTOR_H
