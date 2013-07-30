#include "combiobjectand.h"

CombiObjectAnd::CombiObjectAnd(ObjectDefinition* object_a, ObjectDefinition* object_b)
  : CombiObject(object_a, object_b)
{
}


bool CombiObjectAnd::evalBool()
{
  bool a_inside = m_ObjectA->evalBool();
  bool b_inside = m_ObjectB->evalBool();

  return a_inside && b_inside;
}
