#ifndef CYLINDEROBJECT_H
#define CYLINDEROBJECT_H

class CylinderObject;

#include "objectdefinition.h"

class CylinderObject : public ObjectDefinition
{

protected: // attributes

  /// @todo define wether to use o-coords, or generalized parent coords xyz
  vec3_t m_BottomO; // point on center of bottom (also on axis)
  vec3_t m_AxisO;   // vector from center of bottom to center of top
  real m_Radius;    // radius

  real m_QLength;   // computed
  real m_QRadius;   // computed
  real m_QLR;       // computes L^2 * R^2

public:

    CylinderObject();

    void setParams (real xo_bottom, real yo_bottom, real zo_bottom,
                    real axis_xo, real axis_yo, real axis_zo,
                    real radius);

    virtual bool isInside (const real& xo, const real& yo, const real& zo);

};

#endif // CYLINDEROBJECT_H
