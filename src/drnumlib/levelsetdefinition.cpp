#include "levelsetdefinition.h"

LevelSetDefinition::LevelSetDefinition()
{
}

void LevelSetDefinition::getLowestLevelSets(vector<LevelSetDefinition*>& my_lowest_levelsets)
{
  my_lowest_levelsets.clear();
  my_lowest_levelsets.push_back(this);
}
