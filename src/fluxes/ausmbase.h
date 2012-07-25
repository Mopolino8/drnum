#ifndef AUSMBASE_H
#define AUSMBASE_H

#include "fluxes/compressibleflux.h"

template <class TReconstruction>
class AusmBase : public CompressibleFlux<TReconstruction>
{

protected: // methods

  real M1(real M,real s);
  real M2(real M,real s);
  real M4(real M,real s);
  real P5(real M,real s);

};

template <class TReconstruction>
inline real AusmBase<TReconstruction>::M1(real M,real s)
{
  return 0.5*(M + s*fabs(M));
}

template <class TReconstruction>
inline real AusmBase<TReconstruction>::M2(real M,real s)
{
  return 0.25*s*sqr(M + s);
}

template <class TReconstruction>
inline real AusmBase<TReconstruction>::M4(real M,real s)
{
  if (fabs(M) >= 1) {
    return M1(M, s);
  }
  return M2(M, s)*(1 - 2*s*M2(M, -s));
}

template <class TReconstruction>
inline real AusmBase<TReconstruction>::P5(real M, real s)
{
  if (fabs(M) >= 1) {
    return M1(M,s)/M;
  }
  return M2(M, s)*((2*s - M) - 3*s*M2(M, -s));
}



#endif // AUSMBASE_H
