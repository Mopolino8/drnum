#ifndef BASEFLOWFLUX_H
#define BASEFLOWFLUX_H

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/minmod.h"
#include "reconstruction/secondorder.h"
#include "fluxes/knp.h"
#include "fluxes/vanleer.h"
#include "fluxes/ausmplus.h"
#include "fluxes/ausmdv.h"
//#include "fluxes/ausm.h"
#include "fluxes/roe.h"
#include "fluxes/compressibleslipflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "fluxes/compressiblesmagorinskyflux.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"

class BaseFlowFlux
{

  typedef Roe<5, reconstruction_t, PerfectGas> euler_t;
  typedef CompressibleSmagorinskyFlux<5, 2400, PerfectGas> viscous_t;
  typedef CompressibleSlipFlux<5, reconstruction_t, PerfectGas> wall_t;
  typedef CompressibleFarfieldFlux<5, Upwind1<5>, PerfectGas> farfield_t;

  reconstruction_t m_Reconstruction;
  euler_t          m_EulerFlux;
  viscous_t        m_ViscFlux;
  wall_t           m_WallFlux;
  farfield_t       m_FarfieldFlux;

  real   m_U;
  real   m_P;
  real   m_T;
  real   m_Radius;
  int    m_BasePatch;


public: // methods

  BaseFlowFlux();
  BaseFlowFlux(const BaseFlowFlux& F);

  void setup(real u, real p, real T, real radius, int jet_patch);

  template <typename PATCH> CUDA_DH void xField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

  template <typename PATCH> CUDA_DH void xWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

  template <typename PATCH> CUDA_DH void splitFlux(PATCH *P, splitface_t sf, real* flux);

};


BaseFlowFlux::BaseFlowFlux()
{
}

BaseFlowFlux::BaseFlowFlux(const BaseFlowFlux& F)
{
  m_U = F.m_U;
  m_P = F.m_P;
  m_T = F.m_T;
  m_Radius = F.m_Radius;
  m_BasePatch = F.m_BasePatch;
}

void BaseFlowFlux::setup(real u, real p, real T, real radius, int jet_patch)
{
  m_U = u;
  m_P = p;
  m_T = T;
  m_Radius = radius;
  m_BasePatch = jet_patch;
}


template <typename PATCH>
inline void BaseFlowFlux::xField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.xField(patch, i, j, k, x, y, z, A, flux);
  m_ViscFlux.xField(patch, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void BaseFlowFlux::yField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.yField(patch, i, j, k, x, y, z, A, flux);
  m_ViscFlux.yField(patch, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void BaseFlowFlux::zField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.zField(patch, i, j, k, x, y, z, A, flux);
  m_ViscFlux.zField(patch, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void BaseFlowFlux::xWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.xWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void BaseFlowFlux::yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.yWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void BaseFlowFlux::zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void BaseFlowFlux::xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  dim_t<5> dim;
  real var[5];
  real p, T, u, v, w;
  P->getVar(dim, 0, i, j, k, var);
  PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
  real y0 = 0.5*P->sizeJ()*P->dy();
  real z0 = 0.5*P->sizeK()*P->dz();
  real rr = sqr(y - y0) + sqr(z - z0);
  if (rr < m_Radius*m_Radius && P->getIndex() == m_BasePatch) {
    u = 0;
    v = 0;
    w = 0;
  } else {
    u = m_U;
    v = 0;
    w = 0;
    T = m_T;
    p = m_P;
  }
  flux[0] = A*p/(PerfectGas::R()*T)*u;
  flux[1] = flux[0]*u + A*p;
  flux[2] = flux[0]*v;
  flux[3] = flux[0]*w;
  flux[4] = flux[0]*(PerfectGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

template <typename PATCH>
inline void BaseFlowFlux::yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.yWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void BaseFlowFlux::zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallM(P, i, j, k, x, y, z, A, flux);
}

#endif // BASEFLOWFLUX_H
