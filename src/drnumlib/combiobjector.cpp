#include "combiobjector.h"

CombiObjectOr::CombiObjectOr(ObjectDefinition* object_a, ObjectDefinition* object_b)
  : CombiObject(object_a, object_b)
{
}


CombiObjectOr::CombiObjectOr(ObjectDefinition* object_a)
  : CombiObject(object_a)
{
}


bool CombiObjectOr::evalBool()
{
  bool inside = false;
  for (size_t i_o = 0; i_o < m_Objects.size(); i_o++) {
    if(m_Objects[i_o]->evalBool()) {
      inside = true;
      break;
    }
  }
  return inside;




//  bool a_inside = m_ObjectA->evalBool();
//  bool b_inside = m_ObjectB->evalBool();

//  return a_inside || b_inside;
}
