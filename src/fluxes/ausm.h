#ifndef AUSM_H
#define AUSM_H

#include "fluxes/ausmbase.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class Ausm : public AusmBase
{

  TReconstruction* m_Reconstruction;

public: // methods

  Ausm(TReconstruction* reconstruction) { m_Reconstruction = reconstruction; }

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
void Ausm<TReconstruction, TGas>::xField(CartesianPatch *patch,
                                         size_t i, size_t j, size_t k,
                                         real x, real y, real z,
                                         real A, real* flux)
{
  AUSM_LEFT_PROJX;
  AUSM_RIGHT_PROJX;

  real a   = 0.5*(a_l + a_r);
  real u   = 0.5*(u_l + u_r);
  real M   = u/a;
  real k_l = max(0.0, min(1.0, 0.5*(M + 1)));
  real k_r = 1 - k_l;
  real p   = k_l*p_l + k_r*p_r;
  real Vf  = u*A;
  countFlops(12);

  if (M > 0) {
    flux[0] += Vf*r_l;
    flux[1] += Vf*ru_l;
    flux[2] += Vf*rv_l;
    flux[3] += Vf*rw_l;
    flux[4] += Vf*r_l*H_l;
  } else {
    flux[0] += Vf*r_r;
    flux[1] += Vf*ru_r;
    flux[2] += Vf*rv_r;
    flux[3] += Vf*rw_r;
    flux[4] += Vf*r_r*H_r;
  }
  flux[1] += A*p;
  countFlops(11);
}

template <typename TReconstruction, typename TGas>
void Ausm<TReconstruction, TGas>::yField(CartesianPatch *patch,
                                         size_t i, size_t j, size_t k,
                                         real x, real y, real z,
                                         real A, real* flux)
{
  AUSM_LEFT_PROJY;
  AUSM_RIGHT_PROJY;

  real a   = 0.5*(a_l + a_r);
  real v   = 0.5*(v_l + v_r);
  real M   = v/a;
  real k_l = min(0.0, max(1.0, 0.5*(M + 1)));
  real k_r = 1 - k_l;
  real p   = k_l*p_l + k_r*p_r;
  real Vf  = v*A;
  countFlops(12);

  if (M > 0) {
    flux[0] += Vf*r_l;
    flux[1] += Vf*ru_l;
    flux[2] += Vf*rv_l;
    flux[3] += Vf*rw_l;
    flux[4] += Vf*r_l*H_l;
  } else {
    flux[0] += Vf*r_r;
    flux[1] += Vf*ru_r;
    flux[2] += Vf*rv_r;
    flux[3] += Vf*rw_r;
    flux[4] += Vf*r_r*H_r;
  }
  flux[2] += A*p;
  countFlops(11);
}

template <typename TReconstruction, typename TGas>
void Ausm<TReconstruction, TGas>::zField(CartesianPatch *patch,
                                         size_t i, size_t j, size_t k,
                                         real x, real y, real z,
                                         real A, real* flux)
{
  AUSM_LEFT_PROJZ;
  AUSM_RIGHT_PROJZ;

  real a   = 0.5*(a_l + a_r);
  real w   = 0.5*(w_l + w_r);
  real M   = w/a;
  real k_l = min(0.0, max(1.0, 0.5*(M + 1)));
  real k_r = 1 - k_l;
  real p   = k_l*p_l + k_r*p_r;
  real Vf  = w*A;
  countFlops(12);

  if (M > 0) {
    flux[0] += Vf*r_l;
    flux[1] += Vf*ru_l;
    flux[2] += Vf*rv_l;
    flux[3] += Vf*rw_l;
    flux[4] += Vf*r_l*H_l;
  } else {
    flux[0] += Vf*r_r;
    flux[1] += Vf*ru_r;
    flux[2] += Vf*rv_r;
    flux[3] += Vf*rw_r;
    flux[4] += Vf*r_r*H_r;
  }
  flux[3] += A*p;
  countFlops(11);
}

#endif // AUSM_H
