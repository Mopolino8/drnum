#ifndef CONELEVELSET_H
#define CONELEVELSET_H

class ConeLevelSet;

#include "levelsetdefinition.h"

class ConeLevelSet : public LevelSetDefinition
{

protected: // attributes

  /// @todo define wether to use o-coords, or generalized parent coords xyz
  vec3_t m_BottomO;    // point on center of bottom (also on axis)
  vec3_t m_AxisO;      // vector from center of bottom to center of top
  real m_RadiusBottom; // radius at bottom
  real m_RadiusTop;    // radius at top

  real m_Length;       // computed
  real m_Slope;        // computed
  real m_SlopeCorrect; // computed
  vec3_t m_AxisO_norm; // computed


  real m_QLength;    // computed
  real m_InvQLength; // computed
//  real m_Length;     // computed
  real m_InvLength;  // computed


  real m_QRadius;   // computed
  real m_QLR;       // computes L^2 * R^2


public:

    ConeLevelSet();

    void setParams (real xo_bottom, real yo_bottom, real zo_bottom,
                    real axis_xo, real axis_yo, real axis_zo,
                    real radius_bottom, real radius_top);

    virtual real calcDistance(const real& xo, const real& yo, const real& zo);

};

#endif // CONELEVELSET_H
