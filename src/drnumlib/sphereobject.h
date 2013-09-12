#ifndef SPHEREOBJECT_H
#define SPHEREOBJECT_H

class SphereObject;

#include "objectdefinition.h"

class SphereObject : public ObjectDefinition
{

protected:

  /// @todo define wether to use o-coords, or generalized parent coords xyz
  real m_XoCenter;
  real m_YoCenter;
  real m_ZoCenter;
  real m_Radius;
  real m_QRadius;

public:

  SphereObject();

  void setParams (real xo_center, real yo_center, real zo_center, real radius);

  virtual bool isInside (const real& xo, const real& yo, const real& zo);

};

#endif // SPHEREOBJECT_H
