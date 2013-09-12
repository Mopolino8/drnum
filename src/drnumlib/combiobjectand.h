#ifndef COMBIOBJECTAND_H
#define COMBIOBJECTAND_H

class CombiObjectAnd;

#include "combiobject.h"

class CombiObjectAnd : public CombiObject
{
public:
    CombiObjectAnd(ObjectDefinition* object_a, ObjectDefinition* object_b);
    CombiObjectAnd(ObjectDefinition* object_a);
    virtual bool evalBool();
};

#endif // COMBIOBJECTAND_H
