#ifndef BLOCKCFD_H
#define BLOCKCFD_H

#include <algorithm>
#include <cmath>
#include <iostream>

using namespace std;

typedef double real;

#define restrict __restrict__

inline real sqr(const real &x)
{
  return x*x;
}

inline real sign1(const real &x)
{
  return x >= 0;
}

inline real nonZero(const real &x, const real &eps)
{
  return x + eps*(2*sign1(x) - 1.0)*(1.0 - fabs(tanh(100*x)));
}

#define WITH_VTK
#define VTK_USE_ANSI_STDLIB

#define BUG                             \
  cout << "This seems to be a bug!\n";  \
  cout << "  file: " + __FILE__ + "\n"; \
  cout << "  line: " + __LINE__ + "\n"; \
  exit(EXIT_FAILURE)                    \


#endif // BLOCKCFD_H
