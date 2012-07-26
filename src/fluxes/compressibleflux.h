#ifndef COMPRESSIBLEFLUX_H
#define COMPRESSIBLEFLUX_H

#include "compressiblecartesianpatch.h"

template <class TReconstruction>
class CompressibleFlux : public CompressibleObject
{

public:

  void xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);

};


template <class TReconstruction>
inline void CompressibleFlux<TReconstruction>::xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r  = TReconstruction::xp(P, 0, 0, i, j, k);
  real ir = CHECKED_REAL(1.0/r);
  real ru = TReconstruction::xp(P, 0, 1, i, j, k);
  real rv = TReconstruction::xp(P, 0, 2, i, j, k);
  real rw = TReconstruction::xp(P, 0, 3, i, j, k);
  real u  = ru*ir;
  real v  = rv*ir;
  real w  = rw*ir;
  real rE = TReconstruction::xp(P, 0, 4, i, j, k);
  real T  = CHECKED_REAL((rE*ir - 0.5*(u*u + v*v + w*w))/gasCv());
  real p  = r*gasR()*T;
  flux.var[1] += A*p;
  countFlops(17);
}

template <class TReconstruction>
inline void CompressibleFlux<TReconstruction>::xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r  = TReconstruction::xm(P, 0, 0, i, j, k);
  real ir = CHECKED_REAL(1.0/r);
  real ru = TReconstruction::xm(P, 0, 1, i, j, k);
  real rv = TReconstruction::xm(P, 0, 2, i, j, k);
  real rw = TReconstruction::xm(P, 0, 3, i, j, k);
  real u  = ru*ir;
  real v  = rv*ir;
  real w  = rw*ir;
  real rE = TReconstruction::xm(P, 0, 4, i, j, k);
  real T  = CHECKED_REAL((rE*ir - 0.5*(u*u + v*v + w*w))/gasCv());
  real p  = r*gasR()*T;
  flux.var[1] += A*p;
  countFlops(17);
}


template <class TReconstruction>
inline void CompressibleFlux<TReconstruction>::yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r  = TReconstruction::yp(P, 0, 0, i, j, k);
  real ir = CHECKED_REAL(1.0/r);
  real ru = TReconstruction::yp(P, 0, 1, i, j, k);
  real rv = TReconstruction::yp(P, 0, 2, i, j, k);
  real rw = TReconstruction::yp(P, 0, 3, i, j, k);
  real u  = ru*ir;
  real v  = rv*ir;
  real w  = rw*ir;
  real rE = TReconstruction::yp(P, 0, 4, i, j, k);
  real T  = CHECKED_REAL((rE*ir - 0.5*(u*u + v*v + w*w))/gasCv());
  real p  = r*gasR()*T;
  flux.var[2] += A*p;
  countFlops(17);
}

template <class TReconstruction>
inline void CompressibleFlux<TReconstruction>::yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r  = TReconstruction::ym(P, 0, 0, i, j, k);
  real ir = CHECKED_REAL(1.0/r);
  real ru = TReconstruction::ym(P, 0, 1, i, j, k);
  real rv = TReconstruction::ym(P, 0, 2, i, j, k);
  real rw = TReconstruction::ym(P, 0, 3, i, j, k);
  real u  = ru*ir;
  real v  = rv*ir;
  real w  = rw*ir;
  real rE = TReconstruction::ym(P, 0, 4, i, j, k);
  real T  = CHECKED_REAL((rE*ir - 0.5*(u*u + v*v + w*w))/gasCv());
  real p  = r*gasR()*T;
  flux.var[2] += A*p;
  countFlops(17);
}


template <class TReconstruction>
inline void CompressibleFlux<TReconstruction>::zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r  = TReconstruction::zp(P, 0, 0, i, j, k);
  real ir = CHECKED_REAL(1.0/r);
  real ru = TReconstruction::zp(P, 0, 1, i, j, k);
  real rv = TReconstruction::zp(P, 0, 2, i, j, k);
  real rw = TReconstruction::zp(P, 0, 3, i, j, k);
  real u  = ru*ir;
  real v  = rv*ir;
  real w  = rw*ir;
  real rE = TReconstruction::zp(P, 0, 4, i, j, k);
  real T  = CHECKED_REAL((rE*ir - 0.5*(u*u + v*v + w*w))/gasCv());
  real p  = r*gasR()*T;
  flux.var[3] += A*p;
  countFlops(17);
}

template <class TReconstruction>
inline void CompressibleFlux<TReconstruction>::zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r  = TReconstruction::zm(P, 0, 0, i, j, k);
  real ir = CHECKED_REAL(1.0/r);
  real ru = TReconstruction::zm(P, 0, 1, i, j, k);
  real rv = TReconstruction::zm(P, 0, 2, i, j, k);
  real rw = TReconstruction::zm(P, 0, 3, i, j, k);
  real u  = ru*ir;
  real v  = rv*ir;
  real w  = rw*ir;
  real rE = TReconstruction::zm(P, 0, 4, i, j, k);
  real T  = CHECKED_REAL((rE*ir - 0.5*(u*u + v*v + w*w))/gasCv());
  real p  = r*gasR()*T;
  flux.var[3] += A*p;
  countFlops(17);
}

#endif // COMPRESSIBLEFLUX_H
