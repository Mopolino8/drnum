#ifndef CYLINDERY_H
#define CYLINDERY_H

#include "shape.h"

class CylinderY : public Shape
{

private: // attributes

  real m_Xcentre;
  real m_Ycentre;
  real m_Zcentre;
  real m_Radius;
  real m_XcentreOrig;
  real m_YcentreOrig;
  real m_ZcentreOrig;
  real m_RadiusOrig;


public: // methods

  CylinderY();

  void setCentre(real x, real z);
  void setRadius(real r) { m_Radius = r; }

  bool getBoundaryMetric(real x1, real y1, real z1,
                         real x2, real y2, real z2,
                         real &k, real &nx, real &ny, real &nz);
  bool isInside(real x, real, real z);

  virtual void transform(const Transformation &transformation);
  virtual void reset();

};


inline CylinderY::CylinderY()
{
  setCentre(0, 0);
  setRadius(1.0);
}

inline void CylinderY::setCentre(real x, real z)
{
  m_Xcentre = x;
  m_Ycentre = 0;
  m_Zcentre = z;
}

inline bool CylinderY::getBoundaryMetric(real x1, real y1, real z1,
                                         real x2, real y2, real z2,
                                         real &k, real &nx, real &ny, real &nz)
{
  real dr1 = CHECKED_REAL(sqrt(sqr(x1 - m_Xcentre) + sqr(z1 - m_Zcentre)) - m_Radius);
  real dr2 = CHECKED_REAL(sqrt(sqr(x2 - m_Xcentre) + sqr(z2 - m_Zcentre)) - m_Radius);
  countSqrts(2);
  countFlops(8);

  if (dr1*dr2 > -1e-6*m_Radius) {
    return false;
  }
  countFlops(1);

  k = CHECKED_REAL(fabs(dr1)/(fabs(dr1) + fabs(dr2)));
  countFlops(4);

  real x = k*x2 + (1-k)*x1;
  real y = k*y2 + (1-k)*y1;
  real z = k*z2 + (1-k)*z1;
  countFlops(12);

  nx = x - m_Xcentre;
  ny = 0;
  nz = z - m_Zcentre;
  real nabs = CHECKED_REAL(sqrt(nx*nx + nz*nz));
  nx /= nabs;
  nz /= nabs;
  countSqrts(1);
  countFlops(8);

  return true;
}

inline bool CylinderY::isInside(real x, real, real z)
{
  real dr = CHECKED_REAL(sqrt(sqr(x - m_Xcentre) + sqr(z - m_Zcentre)) - m_Radius);
  countSqrts(1);
  countFlops(4);
  return dr < 0;
}

#endif // CYLINDERY_H
