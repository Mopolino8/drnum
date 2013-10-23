#include "sherelevelset.h"


SphereLevelSet::SphereLevelSet() :
  LevelSetDefinition()
{
}


void SphereLevelSet::setParams (real xo_center, real yo_center, real zo_center, real radius)
{
  m_XoCenter = xo_center;
  m_YoCenter = yo_center;
  m_ZoCenter = zo_center;
  m_Radius = radius;
  m_QRadius = radius * radius;
}


real SphereLevelSet::calcDistance (const real &xo, const real &yo, const real &zo)
{
  // all on xyzo-coords

  real dx = xo - m_XoCenter;
  real dy = yo - m_YoCenter;
  real dz = zo - m_ZoCenter;
  real q_radius_xyzo = dx*dx + dy*dy + dz*dz;
  real radius_xyzo = sqrt(q_radius_xyzo);
  real distance = radius_xyzo - m_Radius;

  return distance;
}
