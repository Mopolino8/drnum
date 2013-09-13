#include "combiobjectand.h"

CombiObjectAnd::CombiObjectAnd(ObjectDefinition* object_a, ObjectDefinition* object_b)
  : CombiObject(object_a, object_b)
{
}


CombiObjectAnd::CombiObjectAnd(ObjectDefinition* object_a)
  : CombiObject(object_a)
{
}


bool CombiObjectAnd::evalBool()
{
  bool inside = true;
  for (size_t i_o = 0; i_o < m_Objects.size(); i_o++) {
    if(!m_Objects[i_o]->evalBool()) {
      inside = false;
      break;
    }
  }
  return inside;





//  bool a_inside = m_ObjectA->evalBool();
//  bool b_inside = m_ObjectB->evalBool();

//  return a_inside && b_inside;
}