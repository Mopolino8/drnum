#ifndef AUSMTOOLS_H
#define AUSMTOOLS_H

#include "blockcfd.h"

class AusmTools
{

public: // static methods

  static real M1(real M,real s);
  static real M2(real M,real s);
  static real M4(real M,real s);
  static real P5(real M,real s);

};

inline real AusmTools::M1(real M,real s)
{
  return 0.5*(M+s*fabs(M));
}

inline real AusmTools::M2(real M,real s)
{
  return 0.25*s*sqr(M+s);
}

inline real AusmTools::M4(real M,real s)
{
  real w = sign1(fabs(M) - 1);
  return w*M1(M, s) + (1 - w)*M2(M, s)*(1 - 2*s*M2(M, -s));
}

inline real AusmTools::P5(real M, real s)
{
  real w = sign1(fabs(M) - 1);
  return w*M1(M, s)/M + M2(M, s)*((2*s - M)-3*s*M*M2(M, -s));
}

#endif // AUSMTOOLS_H
