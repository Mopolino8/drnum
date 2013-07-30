#include "combiobject.h"


CombiObject::CombiObject(ObjectDefinition* object_a, ObjectDefinition* object_b)
{
  m_ObjectA = object_a;
  m_ObjectB = object_b;
  findLowestObjects();
}


void CombiObject::concatLowestObjects(vector<ObjectDefinition*>& other_lowest_objects)
{
  // Append elements of other_lowest_objects onto own list, then sort and make unique
  //.. append
  for (size_t i2 = 0; i2 < other_lowest_objects.size(); i2++) {
    m_LowestObjects.push_back(other_lowest_objects[i2]);
  }
  //.. sort
  sort(m_LowestObjects.begin(), m_LowestObjects.end());
  //.. remove duplicates and resize
  vector<ObjectDefinition*>::iterator it;
  it = unique(m_LowestObjects.begin(), m_LowestObjects.end());
  m_LowestObjects.resize(it - m_LowestObjects.begin());
}


void CombiObject::findLowestObjects()
{
  m_LowestObjects.clear();

  //.. Let object A write directly into own list of "this"
  m_ObjectA->getLowestObjects(m_LowestObjects);

  //.. Let object B write into help list and concatenate on m_LowestObjects
  vector<ObjectDefinition*> your_lowest_objects;
  m_ObjectB->getLowestObjects(your_lowest_objects);
  concatLowestObjects(your_lowest_objects);
}


bool CombiObject::isInside(const real& xo, const real& yo, const real& zo)
{
  // Check wether point is inside lowest objects
  for (size_t i_low = 0; i_low < m_LowestObjects.size(); i_low++) {
    bool low_inside = m_LowestObjects[i_low]->isInside(xo, yo, zo);
    //    m_LowObjResponse[i_low] = low_inside;
    m_LowestObjects[i_low]->setKnownInside(low_inside);
  }

  // Own operation
  bool inside = evalBool();

  return inside;
}

