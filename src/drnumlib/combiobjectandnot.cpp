#include "combiobjectandnot.h"

CombiObjectAndNot::CombiObjectAndNot(ObjectDefinition* object_a, ObjectDefinition* object_b)
  : CombiObject(object_a, object_b)
{
}


CombiObjectAndNot::CombiObjectAndNot(ObjectDefinition* object_a)
  : CombiObject(object_a)
{
}


bool CombiObjectAndNot::evalBool()
{
  // Ask for: "the first object inside and not any other".

  bool inside = m_Objects[0]->evalBool();

  if (inside) {
    for (size_t i_o = 1; i_o < m_Objects.size(); i_o++) {
      if(m_Objects[i_o]->evalBool()) {
        inside = false;
        break;
      }
    }
  }

  return inside;


//  bool a_inside = m_ObjectA->evalBool();
//  bool b_inside = m_ObjectB->evalBool();

//  return a_inside && !b_inside;
}
