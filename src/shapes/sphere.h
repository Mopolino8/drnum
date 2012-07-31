#ifndef SPHERE_H
#define SPHERE_H

#include "blockcfd.h"

class Sphere
{

private: // attributes

  real m_Xcentre;
  real m_Ycentre;
  real m_Zcentre;
  real m_Radius;


public: // methods

  Sphere();

  void setCentre(real x, real y, real z);
  void setRadius(real r) { m_Radius = r; }

  bool getBoundaryMetric(real x1, real y1, real z1,
                         real x2, real y2, real z2,
                         real &k, real &nx, real &ny, real &nz);

};


inline Sphere::Sphere()
{
  setCentre(0, 0, 0);
  setRadius(1.0);
}

inline void Sphere::setCentre(real x, real y, real z)
{
  m_Xcentre = x;
  m_Ycentre = y;
  m_Zcentre = z;
}

inline bool Sphere::getBoundaryMetric(real x1, real y1, real z1,
                                      real x2, real y2, real z2,
                                      real &k, real &nx, real &ny, real &nz)
{
  real dr1 = CHECKED_REAL(sqrt(sqr(x1 - m_Xcentre) + sqr(y1 - m_Ycentre) + sqr(z1 - m_Zcentre)) - m_Radius);
  real dr2 = CHECKED_REAL(sqrt(sqr(x2 - m_Xcentre) + sqr(y2 - m_Ycentre) + sqr(z2 - m_Zcentre)) - m_Radius);
  countSqrts(2);
  countFlops(12);

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
  ny = y - m_Ycentre;
  nz = z - m_Zcentre;
  real nabs = CHECKED_REAL(sqrt(nx*nx + ny*ny + nz*nz));
  nx /= nabs;
  ny /= nabs;
  nz /= nabs;
  countSqrts(1);
  countFlops(11);

  return true;
}


#endif // SPHERE_H
