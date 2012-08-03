#ifndef HALFSPACE_H
#define HALFSPACE_H

#include "shapes/shape.h"

class HalfSpace : public Shape
{

private: // attributes

  real m_X;
  real m_Y;
  real m_Z;
  real m_Nx;
  real m_Ny;
  real m_Nz;

  real m_XOrig;
  real m_YOrig;
  real m_ZOrig;
  real m_NxOrig;
  real m_NyOrig;
  real m_NzOrig;


public: // methods

  HalfSpace();

  void setPoint(real x, real y, real z);
  void setNormal(real nx, real ny, real nz);

  bool getBoundaryMetric(real x1, real y1, real z1,
                         real x2, real y2, real z2,
                         real &k, real &nx, real &ny, real &nz);
  bool isInside(real x, real y, real z);

  virtual void transform(const Transformation &transformation);
  virtual void reset();

};


inline HalfSpace::HalfSpace()
{
  setPoint(0, 0, 0);
  setNormal(0, 0, 1);
}

inline void HalfSpace::setPoint(real x, real y, real z)
{
  m_X = x;
  m_Y = y;
  m_Z = z;
}

inline void HalfSpace::setNormal(real nx, real ny, real nz)
{
  real nabs = sqrt(nx*nx + ny*ny + nz*nz);
  m_Nx = CHECKED_REAL(nx/nabs);
  m_Ny = CHECKED_REAL(ny/nabs);
  m_Nz = CHECKED_REAL(nz/nabs);
  countSqrts(1);
  countFlops(8);
}

inline bool HalfSpace::getBoundaryMetric(real x1, real y1, real z1,
                                      real x2, real y2, real z2,
                                      real &k, real &nx, real &ny, real &nz)
{
  real d1 = (x1 - m_X)*m_Nx + (y1 - m_Y)*m_Ny + (z1 - m_Z)*m_Nz;
  real d2 = (x2 - m_X)*m_Nx + (y2 - m_Y)*m_Ny + (z2 - m_Z)*m_Nz;
  countFlops(16);

  if (d1*d2 > -1e-10) {
    return false;
  }
  countFlops(1);

  k = CHECKED_REAL(fabs(d1)/(fabs(d1) + fabs(d2)));
  countFlops(4);

  nx = m_Nx;
  ny = m_Ny;
  nz = m_Nz;

  return true;
}

inline bool HalfSpace::isInside(real x, real y, real z)
{
  real d = (x - m_X)*m_Nx + (y - m_Y)*m_Ny + (z - m_Z)*m_Nz;
  countFlops(8);
  return d < 0;
}


#endif // HALFSPACE_H
