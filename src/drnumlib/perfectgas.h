// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef PERFECTGAS_H
#define PERFECTGAS_H

#include "cartesianpatch.h"

class PerfectGas
{

public: // static methods

  static CUDA_DH real R(real* = NULL)     { return 287; }       ///< @todo find a concept for this
  static CUDA_DH real gamma(real* = NULL) { return 1.4; }       ///< @todo find a concept for this
  static CUDA_DH real cp(real* = NULL)    { return 1004.5; }    ///< @todo find a concept for this
  static CUDA_DH real cv(real* = NULL)    { return 717.5; }     ///< @todo find a concept for this
  static CUDA_DH real mu(real* = NULL)    { return 1.8e-5; }    ///< @todo find a concept for this
  static CUDA_DH real Pr(real* = NULL)    { return 0.7; }       ///< @todo find a concept for this

  static CUDA_DH void primitiveToConservative(real p, real T, real* var);
  static CUDA_DH void primitiveToConservative(real p, real T, real u, real v, real w, real* var);
  static CUDA_DH void primitiveToConservative(real p, real T, vec3_t U, real* var);
  static CUDA_DH void conservativeToPrimitive(real* var, real& p, real& T);
  static CUDA_DH void conservativeToPrimitive(real* var, real& p, real& T, real& u, real& v, real& w);
  static CUDA_DH void conservativeToPrimitive(real* var, real& p, real& T, vec3_t& U);

};


inline void PerfectGas::primitiveToConservative(real p, real T, real *var)
{
  var[0] = p/(R() * T);
  var[1] = 0;
  var[2] = 0;
  var[3] = 0;
  var[4] = p/(gamma() - 1);
  countFlops(4);
}

inline void PerfectGas::primitiveToConservative(real p, real T, real u, real v, real w, real *var)
{
  var[0] = p/(R() * T);
  var[1] = var[0]*u;
  var[2] = var[0]*v;
  var[3] = var[0]*w;
  var[4] = var[0]*(cv()*T + 0.5*(u*u + v*v + w*w));
  countFlops(14);
}

inline void PerfectGas::primitiveToConservative(real p, real T, vec3_t U, real *var)
{
  primitiveToConservative(p, T, U[0], U[1], U[2], var);
}

inline void PerfectGas::conservativeToPrimitive(real *var, real &p, real &T)
{
  real ir = 1.0/var[0];
  real u  = var[1]*ir;
  real v  = var[2]*ir;
  real w  = var[3]*ir;
  T = (var[4]*ir - 0.5*(u*u + v*v + w*w))/cv();
  p = var[0]*R()*T;
  countFlops(15);
}

inline void PerfectGas::conservativeToPrimitive(real *var, real &p, real &T, real &u, real &v, real &w)
{
  real ir = 1.0/var[0];
  u = var[1]*ir;
  v = var[2]*ir;
  w = var[3]*ir;
  T = (var[4]*ir - 0.5*(u*u + v*v + w*w))/cv();
  p = var[0]*R()*T;
  countFlops(15);
}

inline void PerfectGas::conservativeToPrimitive(real *var, real &p, real &T, vec3_t& U)
{
  conservativeToPrimitive(var, p, T, U[0], U[1], U[2]);
}

#endif // PERFECTGAS_H
