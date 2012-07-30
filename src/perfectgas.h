#ifndef PERFECTGAS_H
#define PERFECTGAS_H

#include "cartesianpatch.h"

class PerfectGas
{

public: // static methods

  static real R(real* = NULL)     { return 287; }    ///< @todo find a concept for this
  static real gamma(real* = NULL) { return 1.4; }    ///< @todo find a concept for this
  static real cp(real* = NULL)    { return 1004.5; } ///< @todo find a concept for this
  static real cv(real* = NULL)    { return 717.5; }  ///< @todo find a concept for this

  static void primitiveToConservative(real p, real T, real* var);
  static void primitiveToConservative(real p, real T, real u, real v, real w, real* var);
  static void conservativeToPrimitive(real* var, real& p, real& T);
  static void conservativeToPrimitive(real* var, real& p, real& T, real& u, real& v, real& w);

};

#endif // PERFECTGAS_H
