#ifndef CARTBOXOBJECT_H
#define CARTBOXOBJECT_H

class CartboxObject;

#include "objectdefinition.h"

class CartboxObject : public ObjectDefinition
{

protected:
  real m_Xo_min;
  real m_Xo_max;
  real m_Yo_min;
  real m_Yo_max;
  real m_Zo_min;
  real m_Zo_max;

public:

  CartboxObject();

  void setParams (real xo_min, real xo_max,
                  real yo_min, real yo_max,
                  real zo_min, real zo_max);

  virtual bool isInside (const real& xo, const real& yo, const real& zo);


};

#endif // CARTBOXOBJECT_H
