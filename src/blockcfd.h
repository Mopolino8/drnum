#ifndef BLOCKCFD_H
#define BLOCKCFD_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ctime>

using namespace std;

typedef double real;

#define RESTRICT __restrict__
#define REGREAL register real

inline real sqr(const real x)
{
  return x*x;
}

inline real sign1(const real x)
{
  return 2.0*(x >= 0) - 1.0;
}

inline real nonZero(const real x, const real eps)
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


#ifdef DEBUG
inline real checkedReal(real x, int line, const char *file_name)
{
  if (isnan(x)) {
    cout << "NaN encountered in file \"" << file_name << "\" at line " << line << endl;
    abort();
  }
  if (isinf(x)) {
    cout << "Inf encountered (division by zero?) in file \"" << file_name << "\" at line " << line << endl;
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


#endif // BLOCKCFD_H
