#include "cartboxobject.h"


CartboxObject::CartboxObject(PatchGrid* patch_grid) :
  BlockObject(patch_grid)
{
}


void CartboxObject::setParams (real xo_min, real xo_max,
                               real yo_min, real yo_max,
                               real zo_min, real zo_max)
{
  m_Xo_min = xo_min;
  m_Xo_max = xo_max;
  m_Yo_min = yo_min;
  m_Yo_max = yo_max;
  m_Zo_min = zo_min;
  m_Zo_max = zo_max;
}


bool CartboxObject::isInside (const real &xo, const real &yo, const real &zo)
{
  bool inside =
      (xo > m_Xo_min) &&
      (xo < m_Xo_max) &&
      (yo > m_Yo_min) &&
      (yo < m_Yo_max) &&
      (zo > m_Zo_min) &&
      (zo < m_Zo_max);

  return inside;
}
