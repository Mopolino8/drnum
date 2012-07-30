#ifndef SPHERE_H
#define SPHERE_H

#include "blockcfd.h"

class Sphere
{

private: // attributes

  real m_X0;
  real m_Y0;
  real m_Z0;
  real m_Radius;

public: // methods

  Sphere();

  void setCentre(real x1, real y1, real z);
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

inline bool Sphere::getBoundaryMetric(real x1, real y1, real z1,
                                      real x2, real y2, real z2,
                                      real &k, real &nx, real &ny, real &nz)
{

}


#endif // SPHERE_H
