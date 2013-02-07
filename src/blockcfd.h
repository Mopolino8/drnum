#ifndef BLOCKCFD_H
#define BLOCKCFD_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ctime>
#include <omp.h>

using namespace std;

#include "namespace_mouse.hh"
#include "mouse_types.hh"    /// @todo mouse_types needed for namespace (?)

using namespace fvmouse;


// configuration macros
//
//#define DEBUG
#define SINGLE_PRECISION



#ifdef SINGLE_PRECISION

typedef float real;

#define FR12 real(0.5)
#define FR13 real(0.3333333)
#define FR23 real(0.6666666)

static real global_eps = 1e-10;

inline int posReal2Int(real v)
{
  v -= 0.5 - 1.5e-7;
  v += 12582912;
  ((char*)&v)[2] &= 63;
  ((char*)&v)[3] &= 0;
  return ((int*)&v)[0];
}

#else

typedef double real;

inline int posReal2Int(real v)
{
  v -= 0.5 - 1.5e-11;
  v += 6755399441055744;
  return ((int*)&v)[0];
}

#endif

#include "math/mathvector.h"
#include "math/smallsquarematrix.h"
#include "debug.h"



using namespace std;


struct real3_t { real x, y, z; };
struct size3_t { size_t i, j, k; };

#define RESTRICT __restrict__
#define REGREAL register real

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


#ifdef DEBUG
inline real checkedReal(real x, int line, const char *file_name)
{
  if (isnan(x)) {
    cout << "NaN encountered in file \"" << file_name << "\" at line " << line << endl;
    GlobalDebug::print();
    abort();
  }
  if (isinf(x)) {
    cout << "Inf encountered (division by zero?) in file \"" << file_name << "\" at line " << line << endl;
    GlobalDebug::print();
    abort();
  }
  return x;
}
#else
inline real checkedReal(real x, int, const char*)
{
  return x;
}
#endif

#define CHECKED_REAL(X) checkedReal(X, __LINE__, __FILE__)


extern unsigned long int global_flops;
extern unsigned long int global_flops_x86;
extern time_t            global_start_time;

#ifdef DEBUG
inline void countFlops(int n)
{
  global_flops     += n;
  global_flops_x86 += n;
}
#else
inline void countFlops(int) {}
#endif

#ifdef DEBUG
inline void countSqrts(int n)
{
  global_flops     += n;
  global_flops_x86 += 15*n;
}
#else
inline void countSqrts(int) {}
#endif

#ifdef DEBUG
inline void countExps(int n)
{
  global_flops     += n;
  global_flops_x86 += 20*n;
}
#else
inline void countExps(int) {}
#endif

#ifdef DEBUG
inline void countLogs(int n)
{
  global_flops     += n;
  global_flops_x86 += 20*n;
}
#else
inline void countLogs(int) {}
#endif

inline void fill(real* var, size_t num_vars, real value)
{
  for (size_t i_var = 0; i_var < num_vars; ++i_var) {
    var[i_var] = value;
  }
}


extern void startTiming();
extern void stopTiming();


inline real sqr(const real x)
{
  return x*x;
  countFlops(1);
}

inline real sign1(const real x)
{
  return 2.0*(x >= 0) - 1.0;
}

inline real nonZero(const real x, const real eps)
{
  if (x < 0) {
    return min(-eps, x);
  } else {
    return max(eps, x);
  }
}

#endif // BLOCKCFD_H
