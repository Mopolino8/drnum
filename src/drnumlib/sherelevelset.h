#ifndef SHERELEVELSET_H
#define SHERELEVELSET_H

class ShereLevelSet;

#include "levelsetdefinition.h"

/**
  * Levelset defintion of a sphere object.
  */
class SphereLevelSet : public LevelSetDefinition
{

protected:

  /// @todo define wether to use o-coords, or generalized parent coords xyz
  real m_XoCenter;
  real m_YoCenter;
  real m_ZoCenter;
  real m_Radius;
  real m_QRadius;

public:

  SphereLevelSet();

  void setParams (real xo_center, real yo_center, real zo_center, real radius);

  virtual real calcDistance (const real& xo, const real& yo, const real& zo);

};

#endif // SHERELEVELSET_H
