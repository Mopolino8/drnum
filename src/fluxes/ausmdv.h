#ifndef AUSMDV_H
#define AUSMDV_H

#include "fluxes/ausmbase.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class AusmDV : public AusmBase
{

  TReconstruction* m_Reconstruction;

public: // methods

  AusmDV(TReconstruction* reconstruction) { m_Reconstruction = reconstruction; }

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
void AusmDV<TReconstruction, TGas>::xField(CartesianPatch *patch,
                                           size_t i, size_t j, size_t k,
                                           real x, real y, real z,
                                           real A, real* flux)
{
  AUSM_LEFT_PROJX;
  AUSM_RIGHT_PROJX;

  real a  = 0.5*(a_l + a_r);
  real fl = p_l/r_l;
  real fr = p_r/r_r;
  real wp = 2*fl/(fl+fr);
  real wm = 2*fr/(fl+fr);
  real Mp = wp*M2(u_l/a, 1) + (1-wp)*M1(u_l/a, 1);
  real Mm = wm*M2(u_r/a,-1) + (1-wm)*M1(u_r/a,-1);
  real p  = P5(u_l/a,1)*p_l + P5(u_r/a,-1)*p_r;
  countFlops(25);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += 0.5*flux[0]*(u_l + u_r) + A*p - 0.5*fabs(flux[0])*(u_r - u_l);
  flux[2] += 0.5*flux[0]*(v_l + v_r)       - 0.5*fabs(flux[0])*(v_r - v_l);
  flux[3] += 0.5*flux[0]*(w_l + w_r)       - 0.5*fabs(flux[0])*(w_r - w_l);
  flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

template <typename TReconstruction, typename TGas>
void AusmDV<TReconstruction, TGas>::yField(CartesianPatch *patch,
                                           size_t i, size_t j, size_t k,
                                           real x, real y, real z,
                                           real A, real* flux)
{
  AUSM_LEFT_PROJY;
  AUSM_RIGHT_PROJY;

  real a  = 0.5*(a_l + a_r);
  real fl = p_l/r_l;
  real fr = p_r/r_r;
  real wp = 2*fl/(fl+fr);
  real wm = 2*fr/(fl+fr);
  real Mp = wp*M2(v_l/a, 1) + (1-wp)*M1(v_l/a, 1);
  real Mm = wm*M2(v_r/a,-1) + (1-wm)*M1(v_r/a,-1);
  real p  = P5(v_l/a,1)*p_l + P5(v_r/a,-1)*p_r;
  countFlops(25);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += 0.5*flux[0]*(u_l + u_r)       - 0.5*fabs(flux[0])*(u_r - u_l);
  flux[2] += 0.5*flux[0]*(v_l + v_r) + A*p - 0.5*fabs(flux[0])*(v_r - v_l);
  flux[3] += 0.5*flux[0]*(w_l + w_r)       - 0.5*fabs(flux[0])*(w_r - w_l);
  flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

template <typename TReconstruction, typename TGas>
void AusmDV<TReconstruction, TGas>::zField(CartesianPatch *patch,
                                           size_t i, size_t j, size_t k,
                                           real x, real y, real z,
                                           real A, real* flux)
{
  AUSM_LEFT_PROJZ;
  AUSM_RIGHT_PROJZ;

  real a  = 0.5*(a_l + a_r);
  real fl = p_l/r_l;
  real fr = p_r/r_r;
  real wp = 2*fl/(fl+fr);
  real wm = 2*fr/(fl+fr);
  real Mp = wp*M2(w_l/a, 1) + (1-wp)*M1(w_l/a, 1);
  real Mm = wm*M2(w_r/a,-1) + (1-wm)*M1(w_r/a,-1);
  real p  = P5(w_l/a,1)*p_l + P5(w_r/a,-1)*p_r;
  countFlops(25);

  flux[0] += a*A*(r_l*Mp + r_r*Mm);
  flux[1] += 0.5*flux[0]*(u_l + u_r)       - 0.5*fabs(flux[0])*(u_r - u_l);
  flux[2] += 0.5*flux[0]*(v_l + v_r)       - 0.5*fabs(flux[0])*(v_r - v_l);
  flux[3] += 0.5*flux[0]*(w_l + w_r) + A*p - 0.5*fabs(flux[0])*(w_r - w_l);
  flux[4] += 0.5*flux[0]*(H_l + H_r)       - 0.5*fabs(flux[0])*(H_r - H_l);
  countFlops(36);
}

#endif // AUSMDV_H
