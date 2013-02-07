#ifndef BOX_H
#define BOX_H

#include "shapes/shape.h"

class Box
{

  real m_X1, m_Y1, m_Z1;
  real m_X2, m_Y2, m_Z2;
  Transformation m_Transform;
  Transformation m_InvTransform;

public: // methods

  Box();

  void setGeometry(real x1, real y1, real z1, real x2, real y2, real z2);

  bool getBoundaryMetric(real x1, real y1, real z1,
                         real x2, real y2, real z2,
                         real &k, real &nx, real &ny, real &nz);
  bool isInside(real x, real y, real z);

  virtual void transform(const Transformation &transformation);
  virtual void reset();

};

#endif // BOX_H
