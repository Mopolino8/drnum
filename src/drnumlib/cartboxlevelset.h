#ifndef CARTBOXLEVELSET_H
#define CARTBOXLEVELSET_H

class CartboxLevelSet;

#include "levelsetdefinition.h"

class CartboxLevelSet : public LevelSetDefinition
{

protected:
  real m_Xo_min;
  real m_Xo_max;
  real m_Yo_min;
  real m_Yo_max;
  real m_Zo_min;
  real m_Zo_max;

public:

  CartboxLevelSet();

  void setParams (real xo_min, real xo_max,
                  real yo_min, real yo_max,
                  real zo_min, real zo_max);

  virtual real calcDistance (const real& xo, const real& yo, const real& zo);


};

#endif // CARTBOXLEVELSET_H
