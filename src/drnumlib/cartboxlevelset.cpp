#include "cartboxlevelset.h"


CartboxLevelSet::CartboxLevelSet() :
  LevelSetDefinition()
{
}


void CartboxLevelSet::setParams (real xo_min, real xo_max,
                                 real yo_min, real yo_max,
                                 real zo_min, real zo_max)
{
  m_Xo_min = xo_min;
  m_Xo_max = xo_max;
  m_Yo_min = yo_min;
  m_Yo_max = yo_max;
  m_Zo_min = zo_min;
  m_Zo_max = zo_max;
}


real CartboxLevelSet::calcDistance (const real &xo, const real &yo, const real &zo)
{
  real distance;
  distance = m_Xo_min - xo;                  // distance from m_Xo_min
  distance = max(distance, (xo - m_Xo_max)); // distance from m_Xo_max
  distance = max(distance, (m_Yo_min - yo)); // distance from m_Yo_min
  distance = max(distance, (yo - m_Yo_max)); // distance from m_Yo_max
  distance = max(distance, (m_Zo_min - zo)); // distance from m_Zo_min
  distance = max(distance, (zo - m_Zo_max)); // distance from m_Zo_max

  return distance;
}
