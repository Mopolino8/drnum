#include "cylinderobject.h"


CylinderObject::CylinderObject() :
  ObjectDefinition()
{
}


void CylinderObject::setParams(real xo_bottom, real yo_bottom, real zo_bottom,
                               real axis_xo, real axis_yo, real axis_zo,
                               real radius)
{
  m_BottomO[0] = xo_bottom;
  m_BottomO[1] = yo_bottom;
  m_BottomO[2] = zo_bottom;

  m_AxisO[0] = axis_xo;
  m_AxisO[1] = axis_yo;
  m_AxisO[2] = axis_zo;

  m_QLength = m_AxisO.abs2();

  m_Radius = radius;
  m_QRadius = radius * radius;

  m_QLR = m_QLength * m_QRadius;
}


bool CylinderObject::isInside (const real &xo, const real &yo, const real &zo)
{
  // all on xyzo-coords

  vec3_t point_xyzo;
  point_xyzo[0] = xo;
  point_xyzo[1] = yo;
  point_xyzo[2] = zo;

  // check between bottom and top
  vec3_t bot2point = point_xyzo - m_BottomO;
  real scal = m_AxisO * bot2point;
  if(scal < 0.) {
    return false;
  }
  else if(scal > m_QLength) {
    return false;
  }

  // check in radial limits
  vec3_t cross = m_AxisO.cross(bot2point);
  real q_dist_q_axilen = cross.abs2();  /// @todo this has unnecessary sqrt()
  if (q_dist_q_axilen > m_QLR) {
    return false;
  }

  // survived until here? => inside
  return true;
}
