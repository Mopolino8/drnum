#include "cylindery.h"

void CylinderY::transform(const Transformation &transformation)
{
  // @todo look at generality (e.g. this would fall over for stretching)
  m_XcentreOrig = m_Xcentre;
  m_YcentreOrig = m_Ycentre;
  m_ZcentreOrig = m_Zcentre;
  m_RadiusOrig = m_Radius;
  real xr = m_Xcentre + m_Radius;
  real yr = m_Ycentre;
  real zr = m_Zcentre;
  transformation(m_Xcentre, m_Ycentre, m_Zcentre);
  transformation(xr, yr, zr);
  m_Radius = sqrt(sqr(xr - m_Xcentre) + sqr(yr - m_Ycentre) + sqr(zr - m_Zcentre));
  countSqrts(1);
  countFlops(6);
}

void CylinderY::reset()
{
  m_Xcentre = m_XcentreOrig;
  m_Ycentre = m_YcentreOrig;
  m_Zcentre = m_ZcentreOrig;
  m_Radius = m_RadiusOrig;
}
