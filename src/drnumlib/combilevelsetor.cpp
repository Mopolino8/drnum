#include "combilevelsetor.h"

CombiLevelSetOr::CombiLevelSetOr(LevelSetDefinition* levelset_a, LevelSetDefinition* levelset_b)
  : CombiLevelSet(levelset_a, levelset_b)
{
}


CombiLevelSetOr::CombiLevelSetOr(LevelSetDefinition* levelset_a)
  : CombiLevelSet(levelset_a)
{
}


real CombiLevelSetOr::evalReal()
{
  // find lowest levelset value in all of m_LevelSets
  real distance = m_LevelSets[0]->evalReal();
  for (size_t i_o = 1; i_o < m_LevelSets.size(); i_o++) {
    real dist_h = m_LevelSets[i_o]->evalReal();
    if(distance > dist_h) {
      distance = dist_h;
    }
  }
  return distance;
}
