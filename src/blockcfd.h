#ifndef BLOCKCFD_H
#define BLOCKCFD_H

#include <algorithm>
#include <cmath>
#include <iostream>

#include "math/mathvector.h"
#include "math/smallsquarematrix.h"

using namespace std;

#include "namespace_mouse.hh"
#include "mouse_types.hh"    /// @todo mouse_types needed for namespace (?)
using namespace fvmouse;

typedef double real;

#define restrict __restrict__

inline real sqr(const real &x)
{
  return x*x;
}

inline real sign1(const real &x)
{
  return 2.0*(x >= 0) - 1.0;
}

inline real nonZero(const real &x, const real &eps)
{
  if (fabs(x) < 0) {
    return min(-eps, x);
  } else {
    return max(eps, x);
  }
  //return x + eps*(2*sign1(x) - 1.0)*(1.0 - fabs(tanh(100*x)));
}

#define WITH_VTK
#define VTK_USE_ANSI_STDLIB

/**
 * Issue an error message and stop the program.
 * Additionally this macro will print the line number and the file
 * where the error occurred.
 */
#define BUG {                             \
  cout << "This seems to be a bug!\n";    \
  cout << "  file: " << __FILE__ << "\n"; \
  cout << "  line: " << __LINE__ << "\n"; \
  abort();                                \
}

#ifdef WITH_VTK
  #include <Qt>
  #ifdef QT_DEBUG
    #define DEBUG
  #endif
#endif


inline real checkedReal(real x, int line, const char *file_name)
{
#ifdef DEBUG
  if (isnan(x)) {
    cout << "NaN encountered in file \"" << file_name << "\" at line " << line << endl;
    abort();
  }
  if (isinf(x)) {
    cout << "Inf encountered (division by zero?) in file \"" << file_name << "\" at line " << line << endl;
    abort();
  }
#endif
  return x;
}

#define CHECKED_REAL(X) checkedReal(X, __LINE__, __FILE__)


extern unsigned long int global_flops;
extern unsigned long int global_flops_x86;

inline void countFlops(int n)
{
#ifdef DEBUG
  global_flops     += n;
  global_flops_x86 += n;
#endif
}

inline void countSqrts(int n)
{
#ifdef DEBUG
  global_flops     += n;
  global_flops_x86 += 15*n;
#endif
}

inline void countExps(int n)
{
#ifdef DEBUG
  global_flops     += n;
  global_flops_x86 += 20*n;
#endif
}

inline void countLogs(int n)
{
#ifdef DEBUG
  global_flops     += n;
  global_flops_x86 += 20*n;
#endif
}


template <unsigned int DIM>
struct RealVec
{
  real var[DIM];
};


#endif // BLOCKCFD_H
