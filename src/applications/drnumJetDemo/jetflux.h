#ifndef JETFLUX_H
#define JETFLUX_H

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/minmod.h"
#include "reconstruction/secondorder.h"
#include "fluxes/knp.h"
#include "fluxes/vanleer.h"
#include "fluxes/ausmplus.h"
//#include "fluxes/ausmdv.h"
#include "fluxes/compressibleslipflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"

class JetFlux
{

  typedef VanLeer<5, reconstruction_t, PerfectGas>              euler_t;
  typedef CompressibleViscFlux<5, PerfectGas>                   viscous_t;
  typedef CompressibleSlipFlux<5, reconstruction_t, PerfectGas> wall_t;
  typedef CompressibleFarfieldFlux<5, Upwind1<5>, PerfectGas>   farfield_t;
  typedef CompressibleFlux<5, PerfectGas>                       split_t;

  reconstruction_t m_Reconstruction;
  euler_t          m_EulerFlux;
  viscous_t        m_ViscFlux;
  wall_t           m_WallFlux;
  farfield_t       m_FarfieldFlux;
  split_t          m_SplitFlux;

  real   m_Ujet;
  real   m_Ufar;
  real   m_Pjet;
  real   m_Pfar;
  real   m_Tjet;
  real   m_Tfar;
  real   m_Radius;
  bool   m_SuperSonic;
  int    m_JetPatch;


public: // methods

  JetFlux();
  JetFlux(const JetFlux& F);

  void setup(real u_jet, real u_far, real p_jet, real p_far, real T_jet, real T_far, real radius, int jet_patch, bool super_sonic);

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


JetFlux::JetFlux()
{
}

JetFlux::JetFlux(const JetFlux& F)
{
  m_Ujet = F.m_Ujet;
  m_Ufar = F.m_Ufar;
  m_Pjet = F.m_Pjet;
  m_Pfar = F.m_Pfar;
  m_Tfar = F.m_Tfar;
  m_Tjet = F.m_Tjet;
  m_Radius = F.m_Radius;
  m_SuperSonic = F.m_SuperSonic;
  m_JetPatch = F.m_JetPatch;
}

void JetFlux::setup(real u_jet, real u_far, real p_jet, real p_far, real T_jet, real T_far, real radius, int jet_patch, bool super_sonic)
{
  m_Ujet = u_jet;
  m_Ufar = u_far;
  m_Pfar = p_far;
  m_Tfar = T_far;
  m_Pjet = p_jet;
  m_Tjet = T_jet;
  m_Radius = radius;
  m_SuperSonic = super_sonic;
  m_JetPatch = jet_patch;
}


template <typename PATCH>
inline void JetFlux::xField
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
inline void JetFlux::yField
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
inline void JetFlux::zField
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
inline void JetFlux::xWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux.xWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.yWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  if (P->getIndex() != m_JetPatch) {
    m_WallFlux.xWallM(P, i, j, k, x, y, z, A, flux);
  } else {
    dim_t<5> dim;
    real var[5];
    real p, T, u, v, w;
    P->getVar(dim, 0, i, j, k, var);
    PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
    real y0 = 0.5*P->sizeJ()*P->dy();
    real z0 = 0.5*P->sizeK()*P->dz();
    real rr = sqr(y - y0) + sqr(z - z0);
    if (rr < m_Radius*m_Radius) {
      u = m_Ujet;
      v = 0;
      w = 0;
      if (m_SuperSonic) {
        p = m_Pjet;
      }
      T = m_Tjet;
    } else {
      u = m_Ufar;
      v = 0;
      w = 0;
      T = m_Tfar;
    }
    flux[0] = A*p/(PerfectGas::R()*T)*u;
    flux[1] = flux[0]*u + A*p;
    flux[2] = flux[0]*v;
    flux[3] = flux[0]*w;
    flux[4] = flux[0]*(PerfectGas::cp()*T + 0.5*(u*u + v*v + w*w));
  }
}

template <typename PATCH>
inline void JetFlux::yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.yWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::splitFlux(PATCH *P, splitface_t sf, real* flux)
{
  m_SplitFlux.splitFlux(P, sf, flux);
}


#endif // JETFLUX_H
