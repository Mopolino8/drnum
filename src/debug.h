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

void print();

#ifdef DEBUG
inline CUDA_DH void xyz(real X, real Y, real Z)
{
  x = X;
  y = Y;
  z = Z;
}
#else
inline void xyz(real, real, real) {}
#endif

#ifdef DEBUG
inline CUDA_DH void ijk(size_t I, size_t J, size_t K)
{
  i = I;
  j = J;
  k = K;
}
#else
inline CUDA_DH void ijk(size_t, size_t, size_t) {}
#endif

}

#endif // DEBUG_H
