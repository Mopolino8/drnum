#include "cylinderlevelset.h"


CylinderLevelSet::CylinderLevelSet() :
  LevelSetDefinition()
{
}


void CylinderLevelSet::setParams(real xo_bottom, real yo_bottom, real zo_bottom,
                                 real axis_xo, real axis_yo, real axis_zo,
                                 real radius)
{
  m_BottomO[0] = xo_bottom;
  m_BottomO[1] = yo_bottom;
  m_BottomO[2] = zo_bottom;

  m_AxisO[0] = axis_xo;
  m_AxisO[1] = axis_yo;
  m_AxisO[2] = axis_zo;

  m_Radius = radius;

  m_Length= m_AxisO.abs();
  m_AxisO_norm = (1./m_Length) * m_AxisO;
}


real CylinderLevelSet::calcDistance (const real &xo, const real &yo, const real &zo)
{
  // all on xyzo-coords

  vec3_t point_xyzo;
  point_xyzo[0] = xo;
  point_xyzo[1] = yo;
  point_xyzo[2] = zo;

  // Compute distance from bottom and top
  vec3_t bot2point = point_xyzo - m_BottomO;
  real scal = m_AxisO_norm * bot2point;
  real dist2bottom = -scal;        // distance from bottom
  real dist2top = scal - m_Length; // distance from top

  // Compute radial distance
  vec3_t cross = m_AxisO_norm.cross(bot2point);
  real dist2surf = cross.abs() - m_Radius;

  // Compute ruling distance
  real distance = dist2bottom;
  if(distance < dist2top) distance = dist2top;
  if(distance < dist2surf) distance = dist2surf;

  return distance;
}
