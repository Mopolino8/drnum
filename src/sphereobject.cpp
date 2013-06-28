#include "sphereobject.h"


SphereObject::SphereObject()
{
}


void SphereObject::setParams (real xo_center, real yo_center, real zo_center, real radius)
{
  m_XoCenter = xo_center;
  m_YoCenter = yo_center;
  m_ZoCenter = zo_center;
  m_Radius = radius;
  m_QRadius = radius * radius;
}


bool SphereObject::isInside (const real &xo, const real &yo, const real &zo)
{
  // all on xyzo-coords

  real dx = xo - m_XoCenter;
  real dy = yo - m_YoCenter;
  real dz = zo - m_ZoCenter;
  real q_dist = dx*dx + dy*dy + dz*dz;

  return (q_dist < m_QRadius);
}
