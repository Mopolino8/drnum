#ifndef AUSMPLUS_H
#define AUSMPLUS_H

#include "fluxes/ausmbase.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class AusmPlus : public AusmBase
{

public: // methods

  void x(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux);
  void y(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux);
  void z(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux);

};


template <typename TReconstruction, typename TGas>
inline void AusmPlus<TReconstruction, TGas>::x(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  AUSM_LEFT_PROJX;
  AUSM_RIGHT_PROJX;

  real a    = 0.5*(a_l + a_r);
  real M    = M4(u_l/a, 1) + M4(u_r/a, -1);
  real Mp   = 0.5*(M + fabs(M));
  real Mm   = 0.5*(M - fabs(M));
  real p    = P5(u_l/a, 1)*p_l + P5(u_r/a, -1)*p_r;
  countFlops(14);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += 0.5*flux[0]*(u_l + u_r) + A*p - 0.5*fabs(flux[0])*(u_r - u_l);
  flux[2] += 0.5*flux[0]*(v_l + v_r)       - 0.5*fabs(flux[0])*(v_r - v_l);
  flux[3] += 0.5*flux[0]*(w_l + w_r)       - 0.5*fabs(flux[0])*(w_r - w_l);
  flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

template <typename TReconstruction, typename TGas>
void AusmPlus<TReconstruction, TGas>::y(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  AUSM_LEFT_PROJY;
  AUSM_RIGHT_PROJY;

  real a    = 0.5*(a_l + a_r);
  real M    = M4(v_l/a, 1) + M4(v_r/a, -1);
  real Mp   = 0.5*(M + fabs(M));
  real Mm   = 0.5*(M - fabs(M));
  real p    = P5(v_l/a, 1)*p_l + P5(v_r/a, -1)*p_r;
  countFlops(14);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += 0.5*flux[0]*(u_l + u_r)       - 0.5*fabs(flux[0])*(u_r - u_l);
  flux[2] += 0.5*flux[0]*(v_l + v_r) + A*p - 0.5*fabs(flux[0])*(v_r - v_l);
  flux[3] += 0.5*flux[0]*(w_l + w_r)       - 0.5*fabs(flux[0])*(w_r - w_l);
  flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

template <typename TReconstruction, typename TGas>
void AusmPlus<TReconstruction, TGas>::z(CartesianPatch *patch, size_t i, size_t j, size_t k, real A, real *flux)
{
  AUSM_LEFT_PROJZ;
  AUSM_RIGHT_PROJZ;

  real a    = 0.5*(a_l + a_r);
  real M    = M4(w_l/a, 1) + M4(w_r/a, -1);
  real Mp   = 0.5*(M + fabs(M));
  real Mm   = 0.5*(M - fabs(M));
  real p    = P5(w_l/a, 1)*p_l + P5(w_r/a, -1)*p_r;
  countFlops(14);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += 0.5*flux[0]*(u_l + u_r)       - 0.5*fabs(flux[0])*(u_r - u_l);
  flux[2] += 0.5*flux[0]*(v_l + v_r)       - 0.5*fabs(flux[0])*(v_r - v_l);
  flux[3] += 0.5*flux[0]*(w_l + w_r) + A*p - 0.5*fabs(flux[0])*(w_r - w_l);
  flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

#endif // AUSMPLUS_H
