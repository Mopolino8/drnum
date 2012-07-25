#include "compressiblecartesianpatch.h"

CompressibleCartesianPatch::CompressibleCartesianPatch()
{
  setNumberOfFields(2);
  setNumberOfVariables(5);
}

void CompressibleCartesianPatch::setState(size_t i, size_t j, size_t k, real p, real T)
{
  f(0, 0, i, j, k) = p/(gasR() * T);
  f(0, 1, i, j, k) = 0;
  f(0, 2, i, j, k) = 0;
  f(0, 3, i, j, k) = 0;
  f(0, 4, i, j, k) = p/(gasGamma() - 1);
}

void CompressibleCartesianPatch::getState(size_t i, size_t j, size_t k, real &p, real &T)
{
  real r  = f(0, 0, i, j, k);
  real ir = 1.0/r;
  real ru = f(0, 1, i, j, k);
  real rv = f(0, 2, i, j, k);
  real rw = f(0, 3, i, j, k);
  real u  = ru*ir;
  real v  = rv*ir;
  real w  = rw*ir;
  real rE = f(0, 4, i, j, k);

  T  = (rE*ir - 0.5*(u*u + v*v + w*w))/gasCv();
  p  = r*gasR()*T;
}

void CompressibleCartesianPatch::getState(size_t i, size_t j, size_t k, real &p, real &u, real &v, real &w, real &T)
{
  real r  = f(0, 0, i, j, k);
  real ir = 1.0/r;
  real ru = f(0, 1, i, j, k);
  real rv = f(0, 2, i, j, k);
  real rw = f(0, 3, i, j, k);
  real rE = f(0, 4, i, j, k);

  u = ru*ir;
  v = rv*ir;
  w = rw*ir;
  T = (rE*ir - 0.5*(u*u + v*v + w*w))/gasCv();
  p = r*gasR()*T;
}

