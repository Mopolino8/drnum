#include "compressiblecartesianpatch.h"

CompressibleCartesianPatch::CompressibleCartesianPatch()
{
  setNumberOfFields(3);
  setNumberOfVariables(5);
}

void CompressibleCartesianPatch::setState(size_t i, size_t j, size_t k, real p, real T)
{
  COMPR_NEW_VARIABLES;

  f(r,  i, j, k) = p/(gasR() * T);
  f(ru, i, j, k) = 0;
  f(rv, i, j, k) = 0;
  f(rw, i, j, k) = 0;
  f(rE, i, j, k) = p/(gasGamma() - 1);
}

void CompressibleCartesianPatch::subStep(real dt)
{
  setFieldToZero(i_res);
  preStep();
  addField(i_old, dt/(dx()*dy()*dz()), i_res, i_new);
  postStep();
}
