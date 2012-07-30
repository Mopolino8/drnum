#ifndef COMPRESSIBLEFARFIELDFLUX_H
#define COMPRESSIBLEFARFIELDFLUX_H

#include "compressiblewallflux.h"

template <typename TReconstruction, typename TGas>
class CompressibleFarfieldFlux
{

  real m_Var[5];

public:

  CompressibleFarfieldFlux() { TGas::primitiveToConservative(1e5, 300, m_Var); }
  void setFarfield(real p, real T, real u, real v, real w) { TGas::primitiveToConservative(p, T, u, v, w, m_Var); }

  void xWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real* flux);
  void xWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real* flux);
  void yWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real* flux);
  void yWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real* flux);
  void zWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real* flux);
  void zWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real* flux);

};


#define FARFIELD_VARS \
  real p, T, u, v, w; \
  real p0, T0, u0, v0, w0; \
  real p1, T1, u1, v1, w1; \
  TGas::conservativeToPrimitive(m_Var, p0, T0, u0, v0, w0); \
  TGas::conservativeToPrimitive(var1, p1, T1, u1, v1, w1); \
  real a0 = CHECKED_REAL(sqrt(TGas::gamma()*TGas::R()*T0)); \
  real a1 = CHECKED_REAL(sqrt(TGas::gamma()*TGas::R()*T1));

#define FARFIELD_SPLIT \
  real M = 0.5*(M0 + M1); \
  if (M >= 0) { \
    if (M >= 1) { \
      p = p0; \
      T = T0; \
      u = u0; \
      v = v0; \
      w = w0; \
    } else { \
      p = p1; \
      T = T0; \
      u = u0; \
      v = v0; \
      w = w0; \
    } \
  } else { \
    if (M <= 1) { \
      p = p1; \
      T = T1; \
      u = u1; \
      v = v1; \
      w = w1; \
    } else { \
      p = p0; \
      T = T1; \
      u = u1; \
      v = v1; \
      w = w1; \
    } \
  }


template <typename TReconstruction, typename TGas>
void CompressibleFarfieldFlux<TReconstruction, TGas>::xWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  real var1[5];
  TReconstruction::project(patch, var1, 0, 5, i-1, j, k, i, j, k);
  FARFIELD_VARS
  real M0 = -CHECKED_REAL(u0/a0);
  real M1 = -CHECKED_REAL(u1/a1);
  FARFIELD_SPLIT
  flux[0] = A*p/(TGas::R()*T)*u;
  flux[1] = flux[0]*u + A*p;
  flux[2] = flux[0]*v;
  flux[3] = flux[0]*w;
  flux[4] = flux[0]*(TGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

template <typename TReconstruction, typename TGas>
void CompressibleFarfieldFlux<TReconstruction, TGas>::xWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  real var1[5];
  TReconstruction::project(patch, var1, 0, 5, i, j, k, i-1, j, k);
  FARFIELD_VARS
  real M0 = CHECKED_REAL(u0/a0);
  real M1 = CHECKED_REAL(u1/a1);
  FARFIELD_SPLIT
  flux[0] = A*p/(TGas::R()*T)*u;
  flux[1] = flux[0]*u + A*p;
  flux[2] = flux[0]*v;
  flux[3] = flux[0]*w;
  flux[4] = flux[0]*(TGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

template <typename TReconstruction, typename TGas>
void CompressibleFarfieldFlux<TReconstruction, TGas>::yWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  real var1[5];
  TReconstruction::project(patch, var1, 0, 5, i, j-1, k, i, j, k);
  FARFIELD_VARS
  real M0 = -CHECKED_REAL(v0/a0);
  real M1 = -CHECKED_REAL(v1/a1);
  FARFIELD_SPLIT
  flux[0] = A*p/(TGas::R()*T)*v;
  flux[1] = flux[0]*u;
  flux[2] = flux[0]*v + A*p;
  flux[3] = flux[0]*w;
  flux[4] = flux[0]*(TGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

template <typename TReconstruction, typename TGas>
void CompressibleFarfieldFlux<TReconstruction, TGas>::yWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  real var1[5];
  TReconstruction::project(patch, var1, 0, 5, i, j, k, i, j-1, k);
  FARFIELD_VARS
  real M0 = CHECKED_REAL(v0/a0);
  real M1 = CHECKED_REAL(v1/a1);
  FARFIELD_SPLIT
  flux[0] = A*p/(TGas::R()*T)*v;
  flux[1] = flux[0]*u;
  flux[2] = flux[0]*v + A*p;
  flux[3] = flux[0]*w;
  flux[4] = flux[0]*(TGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

template <typename TReconstruction, typename TGas>
void CompressibleFarfieldFlux<TReconstruction, TGas>::zWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  real var1[5];
  TReconstruction::project(patch, var1, 0, 5, i, j, k-1, i, j, k);
  FARFIELD_VARS
  real M0 = -CHECKED_REAL(w0/a0);
  real M1 = -CHECKED_REAL(w1/a1);
  FARFIELD_SPLIT
  flux[0] = A*p/(TGas::R()*T)*w;
  flux[1] = flux[0]*u;
  flux[2] = flux[0]*v;
  flux[3] = flux[0]*w + A*p;
  flux[4] = flux[0]*(TGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

template <typename TReconstruction, typename TGas>
void CompressibleFarfieldFlux<TReconstruction, TGas>::zWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  real var1[5];
  TReconstruction::project(patch, var1, 0, 5, i, j, k, i, j, k-1);
  FARFIELD_VARS
  real M0 = CHECKED_REAL(w0/a0);
  real M1 = CHECKED_REAL(w1/a1);
  FARFIELD_SPLIT
  flux[0] = A*p/(TGas::R()*T)*w;
  flux[1] = flux[0]*u;
  flux[2] = flux[0]*v;
  flux[3] = flux[0]*w + A*p;
  flux[4] = flux[0]*(TGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

#endif // COMPRESSIBLEFARFIELDFLUX_H
