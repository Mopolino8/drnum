#ifndef AUSMPLUS_H
#define AUSMPLUS_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class AusmPlus : public CompressibleFlux
{

  TReconstruction* m_Reconstruction;

public: // methods

  AusmPlus(TReconstruction* reconstruction) { m_Reconstruction = reconstruction; }

  void xField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);
  void yField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);
  void zField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);


};


template <typename TReconstruction, typename TGas>
void AusmPlus<TReconstruction, TGas>::xField(CartesianPatch *patch,
                                             size_t i, size_t j, size_t k,
                                             real x, real y, real z,
                                             real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJX;
  COMPRESSIBLE_RIGHT_PROJX;

  real a    = FR12*(a_l + a_r);
  real M    = M4(u_l/a, 1) + M4(u_r/a, -1);
  real Mp   = FR12*(M + fabs(M));
  real Mm   = FR12*(M - fabs(M));
  real p    = P5(u_l/a, 1)*p_l + P5(u_r/a, -1)*p_r;
  countFlops(14);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += FR12*flux[0]*(u_l + u_r) + A*p - FR12*fabs(flux[0])*(u_r - u_l);
  flux[2] += FR12*flux[0]*(v_l + v_r)       - FR12*fabs(flux[0])*(v_r - v_l);
  flux[3] += FR12*flux[0]*(w_l + w_r)       - FR12*fabs(flux[0])*(w_r - w_l);
  flux[4] += FR12*flux[0]*(H_l + H_r)       - FR12*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

template <typename TReconstruction, typename TGas>
void AusmPlus<TReconstruction, TGas>::yField(CartesianPatch *patch,
                                             size_t i, size_t j, size_t k,
                                             real x, real y, real z,
                                             real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJY;
  COMPRESSIBLE_RIGHT_PROJY;

  real a    = FR12*(a_l + a_r);
  real M    = M4(v_l/a, 1) + M4(v_r/a, -1);
  real Mp   = FR12*(M + fabs(M));
  real Mm   = FR12*(M - fabs(M));
  real p    = P5(v_l/a, 1)*p_l + P5(v_r/a, -1)*p_r;
  countFlops(14);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += FR12*flux[0]*(u_l + u_r)       - FR12*fabs(flux[0])*(u_r - u_l);
  flux[2] += FR12*flux[0]*(v_l + v_r) + A*p - FR12*fabs(flux[0])*(v_r - v_l);
  flux[3] += FR12*flux[0]*(w_l + w_r)       - FR12*fabs(flux[0])*(w_r - w_l);
  flux[4] += FR12*flux[0]*(H_l + H_r)       - FR12*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

template <typename TReconstruction, typename TGas>
void AusmPlus<TReconstruction, TGas>::zField(CartesianPatch *patch,
                                             size_t i, size_t j, size_t k,
                                             real x, real y, real z,
                                             real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJZ;
  COMPRESSIBLE_RIGHT_PROJZ;

  real a    = FR12*(a_l + a_r);
  real M    = M4(w_l/a, 1) + M4(w_r/a, -1);
  real Mp   = FR12*(M + fabs(M));
  real Mm   = FR12*(M - fabs(M));
  real p    = P5(w_l/a, 1)*p_l + P5(w_r/a, -1)*p_r;
  countFlops(14);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += FR12*flux[0]*(u_l + u_r)       - FR12*fabs(flux[0])*(u_r - u_l);
  flux[2] += FR12*flux[0]*(v_l + v_r)       - FR12*fabs(flux[0])*(v_r - v_l);
  flux[3] += FR12*flux[0]*(w_l + w_r) + A*p - FR12*fabs(flux[0])*(w_r - w_l);
  flux[4] += FR12*flux[0]*(H_l + H_r)       - FR12*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

#endif // AUSMPLUS_H
