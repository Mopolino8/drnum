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
#ifndef BLOCKCFD_H
#define BLOCKCFD_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <ctime>
#include <cassert>
#include <cstdio>
#include <omp.h>
#include <limits>

using namespace std;


// configuration macros
//
//#define DEBUG
#define SINGLE_PRECISION


#ifdef __CUDACC__
  #define CUDA_DO __device__
  #define CUDA_HO __host__
  #define CUDA_DH __device__ __host__
#else
  #define CUDA_DO
  #define CUDA_HO
  #define CUDA_DH
#endif


#define WITH_VTK
#ifdef SINGLE_PRECISION

typedef float real;

#define MAX_REAL numeric_limits<real>::max()
#define MIN_REAL numeric_limits<real>::min()

#define FR12 real(0.5)
#define FR13 real(0.3333333)
#define FR14 real(0.25)
#define FR23 real(0.6666666)

#define GLOBAL_EPS real(1e-10)

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



#define RESTRICT __restrict__
#define REGREAL real
//register real

#define VTK_USE_ANSI_STDLIB

/**
 * Issue an error message and stop the program.
 * Additionally this macro will print the line number and the file
 * where the error occurred.
 */
#define BUG {                             \
  printf("This seems to be a bug!\n");    \
  printf("  file: %s\n", __FILE__);       \
  printf("  line: %d\n", __LINE__);       \
  assert(false);                          \
  exit(EXIT_FAILURE);                     \
}

/**
 * Issue an error message and stop the program.
 * Additionally this macro will print the line number and the file
 * where the error occurred.
 */
#define ERROR(msg) {                      \
  printf("An error occured: %s\n", msg);  \
  printf("  file: %s\n", __FILE__);       \
  printf("  line: %d\n", __LINE__);       \
  assert(false);                          \
  exit(EXIT_FAILURE);                     \
}

#ifdef WITH_VTK
  #include <QtCore/Qt>
  #ifdef QT_DEBUG
    #ifndef DEBUG
      #define DEBUG
    #endif
  #endif
#endif

#ifdef __CUDACC__
  #ifdef DEBUG
    #undef DEBUG
  #endif
#endif

#ifdef DEBUG
  inline CUDA_DH real checkedReal(real x, int line, const char *file_name)
  {
    if (isnan(x)) {
      //cout << "NaN encountered in file \"" << file_name << "\" at line " << line << endl;
      printf("NaN encountered in file \"%s\" at line %d\n", file_name, line);
      GlobalDebug::print();
      assert(false);
    }
    if (isinf(x)) {
      //cout << "Inf encountered (division by zero?) in file \"" << file_name << "\" at line " << line << endl;
      printf("Inf encountered (division by zero?) in file \"%s\" at line %d\n", file_name, line);
      GlobalDebug::print();
      assert(false);
    }
    return x;
  }
#else
  inline CUDA_DH real checkedReal(real x, int, const char*)
  {
    return x;
  }
#endif

#define CHECKED_REAL(X) checkedReal(X, __LINE__, __FILE__)


extern unsigned long int global_flops;
extern unsigned long int global_flops_x86;
extern time_t            global_start_time;

#ifdef DEBUG
  #ifndef __CUDACC__
    inline void countFlops(int n)
    {
      global_flops     += n;
      global_flops_x86 += n;
    }
  #else
    inline CUDA_DH void countFlops(int) {}
  #endif
#else
  inline CUDA_DH void countFlops(int) {}
#endif

#ifdef DEBUG
  #ifndef __CUDACC__
    inline void countSqrts(int n)
    {
      global_flops     += n;
      global_flops_x86 += 15*n;
    }
  #else
    inline CUDA_DH void countSqrts(int) {}
  #endif
#else
  inline CUDA_DH void countSqrts(int) {}
#endif

#ifdef DEBUG
  #ifndef __CUDACC__
    inline void countExps(int n)
    {
      global_flops     += n;
      global_flops_x86 += 20*n;
    }
  #else
    inline CUDA_DH void countExps(int) {}
  #endif
#else
  inline CUDA_DH void countExps(int) {}
#endif

#ifdef DEBUG
  #ifndef __CUDACC__
    inline void countLogs(int n)
    {
      global_flops     += n;
      global_flops_x86 += 20*n;
    }
  #else
    inline CUDA_DH void countLogs(int) {}
  #endif
#else
  inline CUDA_DH void countLogs(int) {}
#endif

inline CUDA_DH void fill(real* var, size_t num_vars, real value)
{
  for (size_t i_var = 0; i_var < num_vars; ++i_var) {
    var[i_var] = value;
  }
}


extern void startTiming();
extern void stopTiming();
extern void printTiming();


inline CUDA_DH real sqr(const real x)
{
  countFlops(1);
  return x*x;
}

inline CUDA_DH real sign1(const real x)
{
  return 2.0*(x >= 0) - 1.0;
}

inline CUDA_DH real nonZero(const real x, const real eps)
{
  if (x < 0) {
    return min(-eps, x);
  } else {
    return max(eps, x);
  }
}

struct real3_t { real x, y, z; };

struct size3_t
{
  size_t i, j, k;
  CUDA_DH size3_t()                             { i = 0; j = 0; k = 0; }
  CUDA_DH size3_t(size_t i, size_t j, size_t k) { this->i = i; this->j = j; this->k = k; }
  CUDA_DH bool operator==(const size3_t& s)     { return i == s.i && j == s.j && k == s.k; }
};

template <unsigned int DIM>
struct dim_t
{
  static const unsigned int dim = DIM;
  unsigned int operator() () { return DIM; }
};


#endif // BLOCKCFD_H
