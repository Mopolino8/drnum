#ifndef AUSMBASE_H
#define AUSMBASE_H

#include "compressibleobject.h"

class AusmBase : public CompressibleObject
{

protected: // methods

  real M1(real M,real s);
  real M2(real M,real s);
  real M4(real M,real s);
  real P5(real M,real s);

};

inline real AusmBase::M1(real M,real s)
{
  return 0.5*(M + s*fabs(M));
}

inline real AusmBase::M2(real M,real s)
{
  return 0.25*s*sqr(M + s);
}

inline real AusmBase::M4(real M,real s)
{
  if (fabs(M) >= 1) {
    return M1(M, s);
  }
  return M2(M, s)*(1 - 2*s*M2(M, -s));
}

inline real AusmBase::P5(real M, real s)
{
  if (fabs(M) >= 1) {
    return M1(M,s)/M;
  }
  return M2(M, s)*((2*s - M) - 3*s*M2(M, -s));
}

#endif // AUSMBASE_H
