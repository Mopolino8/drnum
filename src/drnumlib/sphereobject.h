#ifndef SPHEREOBJECT_H
#define SPHEREOBJECT_H

class SphereObject;

#include "blockobject.h"

class SphereObject : public BlockObject
{

protected:
  real m_XoCenter;
  real m_YoCenter;
  real m_ZoCenter;
  real m_Radius;
  real m_QRadius;

public:

  SphereObject(PatchGrid* patch_grid);

  void setParams (real xo_center, real yo_center, real zo_center, real radius);

  virtual bool isInside (const real& xo, const real& yo, const real& zo);

};

#endif // SPHEREOBJECT_H
