#include "combilevelsetand.h"

CombiLevelSetAnd::CombiLevelSetAnd(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b)
  : CombiLevelSet(levelset_a, levelset_b)
{
}


CombiLevelSetAnd::CombiLevelSetAnd(LevelSetDefinition* levelset_a)
  : CombiLevelSet(levelset_a)
{
}


real CombiLevelSetAnd::evalReal()
{
  // find highest levelset value in all of m_LevelSets
  real distance = m_LevelSets[0]->evalReal();
  for (size_t i_o = 1; i_o < m_LevelSets.size(); i_o++) {
    real dist_h = m_LevelSets[i_o]->evalReal();
    if(distance < dist_h) {
      distance = dist_h;
    }
  }
  return distance;
}

