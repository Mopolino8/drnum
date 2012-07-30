#include "perfectgas.h"

void PerfectGas::primitiveToConservative(real p, real T, real *var)
{
  var[0] = p/(R() * T);
  var[1] = 0;
  var[2] = 0;
  var[3] = 0;
  var[4] = p/(gamma() - 1);
}

void PerfectGas::primitiveToConservative(real p, real T, real u, real v, real w, real *var)
{
  var[0] = p/(R() * T);
  var[1] = var[0]*u;
  var[2] = var[0]*v;
  var[3] = var[0]*w;
  var[4] = p/(gamma() - 1);
}

void PerfectGas::conservativeToPrimitive(real *var, real &p, real &T)
{
  real ir = 1.0/var[0];
  real u  = var[1]*ir;
  real v  = var[2]*ir;
  real w  = var[3]*ir;
  T = (var[4]*ir - 0.5*(u*u + v*v + w*w))/cv();
  p = var[0]*R()*T;
}

void PerfectGas::conservativeToPrimitive(real *var, real &p, real &T, real &u, real &v, real &w)
{
  real ir = 1.0/var[0];
  u = var[1]*ir;
  v = var[2]*ir;
  w = var[3]*ir;
  T = (var[4]*ir - 0.5*(u*u + v*v + w*w))/cv();
  p = var[0]*R()*T;
}
