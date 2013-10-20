#include "conelevelset.h"

ConeLevelSet::ConeLevelSet() :
  LevelSetDefinition()
{
}


void ConeLevelSet::setParams(real xo_bottom, real yo_bottom, real zo_bottom,
                             real axis_xo, real axis_yo, real axis_zo,
                             real radius_bottom, real radius_top)
{
  m_BottomO[0] = xo_bottom;
  m_BottomO[1] = yo_bottom;
  m_BottomO[2] = zo_bottom;

  m_AxisO[0] = axis_xo;
  m_AxisO[1] = axis_yo;
  m_AxisO[2] = axis_zo;

  m_RadiusBottom = radius_bottom;
  m_RadiusTop = radius_top;
  m_Length = m_AxisO.abs();
  real delta_radius = m_RadiusTop - m_RadiusBottom;
  m_Slope = delta_radius / m_Length;

  real cord_len = sqrt(delta_radius * delta_radius + m_Length * m_Length);
  m_SlopeCorrect = m_Length / cord_len; // = cos(alfa)

  m_AxisO_norm = (1./m_Length) * m_AxisO;
}


real ConeLevelSet::calcDistance (const real &xo, const real &yo, const real &zo)
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
  //.. radial distance of point from axis
  vec3_t cross = m_AxisO_norm.cross(bot2point);
  real dist2axis = cross.abs();
  //.. radial distance of cone surface from axis
  real conerad2axis = m_RadiusBottom + m_Slope * scal;
  real radial_dist = dist2axis - conerad2axis;
  real dist2surf = radial_dist * m_SlopeCorrect;

  // Compute ruling distance
  real distance = dist2bottom;
  if(distance < dist2top) distance = dist2top;
  if(distance < dist2surf) distance = dist2surf;

  return distance;
}
