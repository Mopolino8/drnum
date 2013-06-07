#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>

#include "blockcfd.h"

namespace GlobalDebug
{

using namespace std;

extern size_t i;
extern size_t j;
extern size_t k;
extern real   x;
extern real   y;
extern real   z;

#ifndef __CUDACC__
  inline void CUDA_DH print()
  {
    printf("last variables:\n");
    printf("---------------\n");
    printf("i = %u\n", i);
    printf("j = %u\n", j);
    printf("k = %u\n", k);
    printf("x = %f\n", x);
    printf("y = %f\n", y);
    printf("z = %f\n", z);
    printf("\n");
  }
#else
  inline void CUDA_DH print() {}
#endif


#ifdef DEBUG
  #ifndef __CUDACC__
    inline CUDA_DH void xyz(real X, real Y, real Z)
    {
      x = X;
      y = Y;
      z = Z;
    }
  #else
    inline CUDA_DH void xyz(real, real, real) {}
  #endif
#else
  inline CUDA_DH void xyz(real, real, real) {}
#endif

#ifdef DEBUG
  #ifndef __CUDACC__
    inline CUDA_DH void ijk(size_t I, size_t J, size_t K)
    {
      i = I;
      j = J;
      k = K;
    }
  #else
    inline CUDA_DH void ijk(size_t, size_t, size_t) {}
  #endif
#else
  inline CUDA_DH void ijk(size_t, size_t, size_t) {}
#endif

}

#define DBG_PRT_INT(X) printf("variable %s = %d\n", #X, X);
#define DBG_PRT_REAL(X) printf("variable %s = %f\n", #X, X);

#endif // DEBUG_H
