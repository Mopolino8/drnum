#include "combilevelsetandnot.h"

CombiLevelSetAndNot::CombiLevelSetAndNot(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b)
  : CombiLevelSet(levelset_a, levelset_b)
{
}


CombiLevelSetAndNot::CombiLevelSetAndNot(LevelSetDefinition* levelset_a)
  : CombiLevelSet(levelset_a)
{
}


real CombiLevelSetAndNot::evalReal()
{
  // find highest value of levelset A (stored in m_LevelSets[0]) and
  // negative levelset value of any other (these are the B-objects)
  real distance = m_LevelSets[0]->evalReal();
  for (size_t i_o = 1; i_o < m_LevelSets.size(); i_o++) {
    real dist_h = -m_LevelSets[i_o]->evalReal();
    if(distance < dist_h) {
      distance = dist_h;
    }
  }
  return distance;
}
