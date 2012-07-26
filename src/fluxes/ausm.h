#ifndef AUSMHAENEL_H
#define AUSMHAENEL_H

#include "fluxes/ausmbase.h"
#include "cartesianpatch.h"

template <class TReconstruction>
class Ausm : public AusmBase
{

public: // methods

  void x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);

};


template <class TReconstruction>
inline void Ausm<TReconstruction>::x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r_l  = TReconstruction::xp(P, 0, 0, i, j, k);
  real ir_l = CHECKED_REAL(1.0/r_l);
  real ru_l = TReconstruction::xp(P, 0, 1, i, j, k);
  real rv_l = TReconstruction::xp(P, 0, 2, i, j, k);
  real rw_l = TReconstruction::xp(P, 0, 3, i, j, k);
  real u_l  = ru_l*ir_l;
  real v_l  = rv_l*ir_l;
  real w_l  = rw_l*ir_l;
  real rE_l = TReconstruction::xp(P, 0, 4, i, j, k);
  real T_l  = CHECKED_REAL((rE_l*ir_l - 0.5*(u_l*u_l + v_l*v_l + w_l*w_l))/gasCv());
  real p_l  = r_l*gasR()*T_l;
  real a_l  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_l));
  real H_l  = (rE_l + p_l)/r_l;
  countFlops(19);
  countSqrts(1);

  real r_r  = TReconstruction::xm(P, 0, 0, i, j, k);
  real ir_r = CHECKED_REAL(1.0/r_r);
  real ru_r = TReconstruction::xm(P, 0, 1, i, j, k);
  real rv_r = TReconstruction::xm(P, 0, 2, i, j, k);
  real rw_r = TReconstruction::xm(P, 0, 3, i, j, k);
  real u_r  = ru_r*ir_r;
  real v_r  = rv_r*ir_r;
  real w_r  = rw_r*ir_r;
  real rE_r = TReconstruction::xm(P, 0, 4, i, j, k);
  real T_r  = CHECKED_REAL((rE_r*ir_r - 0.5*(u_r*u_r + v_r*v_r + w_r*w_r))/gasCv());
  real p_r  = r_r*gasR()*T_r;
  real a_r  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_r));
  real H_r  = (rE_r + p_r)/r_r;
  countFlops(19);
  countSqrts(1);

  real a   = 0.5*(a_l + a_r);
  real u   = 0.5*(u_l + u_r);
  real M   = u/a;
  real k_l = max(0.0, min(1.0, 0.5*(M + 1)));
  real k_r = 1 - k_l;
  real p   = k_l*p_l + k_r*p_r;
  real Vf  = u*A;
  countFlops(12);

  if (M > 0) {
    flux.var[0] += Vf*r_l;
    flux.var[1] += Vf*ru_l;
    flux.var[2] += Vf*rv_l;
    flux.var[3] += Vf*rw_l;
    flux.var[4] += Vf*r_l*H_l;
  } else {
    flux.var[0] += Vf*r_r;
    flux.var[1] += Vf*ru_r;
    flux.var[2] += Vf*rv_r;
    flux.var[3] += Vf*rw_r;
    flux.var[4] += Vf*r_r*H_r;
  }
  flux.var[1] += A*p;
  countFlops(11);
}

template <class TReconstruction>
void Ausm<TReconstruction>::y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r_l  = TReconstruction::yp(P, 0, 0, i, j, k);
  real ir_l = CHECKED_REAL(1.0/r_l);
  real ru_l = TReconstruction::yp(P, 0, 1, i, j, k);
  real rv_l = TReconstruction::yp(P, 0, 2, i, j, k);
  real rw_l = TReconstruction::yp(P, 0, 3, i, j, k);
  real u_l  = ru_l*ir_l;
  real v_l  = rv_l*ir_l;
  real w_l  = rw_l*ir_l;
  real rE_l = TReconstruction::yp(P, 0, 4, i, j, k);
  real T_l  = CHECKED_REAL((rE_l*ir_l - 0.5*(u_l*u_l + v_l*v_l + w_l*w_l))/gasCv());
  real p_l  = r_l*gasR()*T_l;
  real a_l  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_l));
  real H_l  = (rE_l + p_l)/r_l;
  countFlops(19);
  countSqrts(1);

  real r_r  = TReconstruction::ym(P, 0, 0, i, j, k);
  real ir_r = CHECKED_REAL(1.0/r_r);
  real ru_r = TReconstruction::ym(P, 0, 1, i, j, k);
  real rv_r = TReconstruction::ym(P, 0, 2, i, j, k);
  real rw_r = TReconstruction::ym(P, 0, 3, i, j, k);
  real u_r  = ru_r*ir_r;
  real v_r  = rv_r*ir_r;
  real w_r  = rw_r*ir_r;
  real rE_r = TReconstruction::ym(P, 0, 4, i, j, k);
  real T_r  = CHECKED_REAL((rE_r*ir_r - 0.5*(u_r*u_r + v_r*v_r + w_r*w_r))/gasCv());
  real p_r  = r_r*gasR()*T_r;
  real a_r  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_r));
  real H_r  = (rE_r + p_r)/r_r;
  countFlops(19);
  countSqrts(1);

  real a   = 0.5*(a_l + a_r);
  real v   = 0.5*(v_l + v_r);
  real M   = v/a;
  real k_l = min(0.0, max(1.0, 0.5*(M + 1)));
  real k_r = 1 - k_l;
  real p   = k_l*p_l + k_r*p_r;
  real Vf  = v*A;
  countFlops(12);

  if (M > 0) {
    flux.var[0] += Vf*r_l;
    flux.var[1] += Vf*ru_l;
    flux.var[2] += Vf*rv_l;
    flux.var[3] += Vf*rw_l;
    flux.var[4] += Vf*r_l*H_l;
  } else {
    flux.var[0] += Vf*r_r;
    flux.var[1] += Vf*ru_r;
    flux.var[2] += Vf*rv_r;
    flux.var[3] += Vf*rw_r;
    flux.var[4] += Vf*r_r*H_r;
  }
  flux.var[2] += A*p;
  countFlops(11);
}

template <class TReconstruction>
void Ausm<TReconstruction>::z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  real r_l  = TReconstruction::zp(P, 0, 0, i, j, k);
  real ir_l = CHECKED_REAL(1.0/r_l);
  real ru_l = TReconstruction::zp(P, 0, 1, i, j, k);
  real rv_l = TReconstruction::zp(P, 0, 2, i, j, k);
  real rw_l = TReconstruction::zp(P, 0, 3, i, j, k);
  real u_l  = ru_l*ir_l;
  real v_l  = rv_l*ir_l;
  real w_l  = rw_l*ir_l;
  real rE_l = TReconstruction::zp(P, 0, 4, i, j, k);
  real T_l  = CHECKED_REAL((rE_l*ir_l - 0.5*(u_l*u_l + v_l*v_l + w_l*w_l))/gasCv());
  real p_l  = r_l*gasR()*T_l;
  real a_l  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_l));
  real H_l  = (rE_l + p_l)/r_l;
  countFlops(19);
  countSqrts(1);

  real r_r  = TReconstruction::zm(P, 0, 0, i, j, k);
  real ir_r = CHECKED_REAL(1.0/r_r);
  real ru_r = TReconstruction::zm(P, 0, 1, i, j, k);
  real rv_r = TReconstruction::zm(P, 0, 2, i, j, k);
  real rw_r = TReconstruction::zm(P, 0, 3, i, j, k);
  real u_r  = ru_r*ir_r;
  real v_r  = rv_r*ir_r;
  real w_r  = rw_r*ir_r;
  real rE_r = TReconstruction::zm(P, 0, 4, i, j, k);
  real T_r  = CHECKED_REAL((rE_r*ir_r - 0.5*(u_r*u_r + v_r*v_r + w_r*w_r))/gasCv());
  real p_r  = r_r*gasR()*T_r;
  real a_r  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_r));
  real H_r  = (rE_r + p_r)/r_r;
  countFlops(19);
  countSqrts(1);

  real a   = 0.5*(a_l + a_r);
  real w   = 0.5*(w_l + w_r);
  real M   = w/a;
  real k_l = min(0.0, max(1.0, 0.5*(M + 1)));
  real k_r = 1 - k_l;
  real p   = k_l*p_l + k_r*p_r;
  real Vf  = w*A;
  countFlops(12);

  if (M > 0) {
    flux.var[0] += Vf*r_l;
    flux.var[1] += Vf*ru_l;
    flux.var[2] += Vf*rv_l;
    flux.var[3] += Vf*rw_l;
    flux.var[4] += Vf*r_l*H_l;
  } else {
    flux.var[0] += Vf*r_r;
    flux.var[1] += Vf*ru_r;
    flux.var[2] += Vf*rv_r;
    flux.var[3] += Vf*rw_r;
    flux.var[4] += Vf*r_r*H_r;
  }
  flux.var[3] += A*p;
  countFlops(11);
}

#endif // AUSMHAENEL_H
