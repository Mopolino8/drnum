#include "objectdefinition.h"

ObjectDefinition::ObjectDefinition()
{
}


void ObjectDefinition::getLowestObjects(vector<ObjectDefinition*>& my_lowest_objects)
{
  my_lowest_objects.clear();
  my_lowest_objects.push_back(this);
}
