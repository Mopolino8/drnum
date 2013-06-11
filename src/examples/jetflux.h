#ifndef JETFLUX_H
#define JETFLUX_H

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/minmod.h"
#include "reconstruction/secondorder.h"
#include "fluxes/knp.h"
//#include "fluxes/vanleer.h"
#include "fluxes/ausmplus.h"
//#include "fluxes/ausmdv.h"
#include "fluxes/compressiblewallflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"

class JetFlux
{

  //typedef Upwind1 reconstruction_t;
  //typedef Upwind2<VanAlbada> reconstruction_t;
  typedef Upwind2<SecondOrder> reconstruction_t;
  typedef AusmPlus<reconstruction_t, PerfectGas> euler_t;
  typedef CompressibleViscFlux<PerfectGas> viscous_t;
  typedef CompressibleWallFlux<reconstruction_t, PerfectGas> wall_t;

  reconstruction_t m_Reconstruction;
  euler_t          m_EulerFlux;
  viscous_t        m_ViscFlux;
  wall_t           m_WallFlux;

  real m_Ujet;
  real m_Ufar;
  real m_Pjet;
  real m_Pfar;
  real m_Tjet;
  real m_Tfar;
  real m_Y1;
  real m_Y2;
  real m_Z1;
  real m_Z2;
  bool m_SuperSonic;


public: // methods

  JetFlux(real u_jet, real u_far, real p_jet, real p_far, real T_jet, real T_far, real y1, real y2, real z1, real z2, bool super_sonic);
  JetFlux();
  JetFlux(const JetFlux& F);

  template <typename PATCH> CUDA_DH void xField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

  template <typename PATCH> CUDA_DH void xWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

};


JetFlux::JetFlux(real u_jet, real u_far, real p_jet, real p_far, real T_jet, real T_far, real y1, real y2, real z1, real z2, bool super_sonic)
{
  m_Ujet = u_jet;
  m_Ufar = u_far;
  m_Pfar = p_far;
  m_Tfar = T_far;
  m_Pjet = p_jet;
  m_Tjet = T_jet;
  m_Y1 = y1;
  m_Y2 = y2;
  m_Z1 = z1;
  m_Z2 = z2;
  m_SuperSonic = super_sonic;
}

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
  m_Y1 = F.m_Y1;
  m_Y2 = F.m_Y2;
  m_Z1 = F.m_Z1;
  m_Z2 = F.m_Z2;
  m_SuperSonic = F.m_SuperSonic;
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
  real var[5];
  real p, T, u, v, w;
  P->getVar(0, i-1, j, k, var);
  PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
  real Ma = u/sqrt(PerfectGas::gamma(var)*PerfectGas::R(var)*T);
  if (Ma < 1) {
    p = m_Pfar;
  }
  flux[0] = A*p/(PerfectGas::R()*T)*u;
  flux[1] = flux[0]*u + A*p;
  flux[2] = flux[0]*v;
  flux[3] = flux[0]*w;
  flux[4] = flux[0]*(PerfectGas::cp()*T + 0.5*(u*u + v*v + w*w));
}

template <typename PATCH>
inline void JetFlux::yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux.yWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux.zWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  /// @todo Oli: modified here to allow shifted patches. Only valid, if aligned to o-coords.
  vec3_t xyzo_ref = P->accessBBoxXYZoMin();
  real yo = y + xyzo_ref[1];
  real zo = z + xyzo_ref[2];

  real var[5];
  real p, T, u, v, w;
  P->getVar(0, i, j, k, var);
  PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
  if (yo > m_Y1 && yo < m_Y2 && zo > m_Z1 && zo < m_Z2) {
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

template <typename PATCH>
inline void JetFlux::yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux.yWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void JetFlux::zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux.zWallM(P, i, j, k, x, y, z, A, flux);
}



#endif // JETFLUX_H
