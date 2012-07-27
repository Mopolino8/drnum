#ifndef AUSMPLUS_H
#define AUSMPLUS_H

#include "fluxes/ausmbase.h"
#include "cartesianpatch.h"

template <class TReconstruction>
class AusmPlus : public AusmBase
{

public: // methods

  void x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);

};


template <class TReconstruction>
inline void AusmPlus<TReconstruction>::x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
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

  real a    = 0.5*(a_l + a_r);
  real M    = AusmTools::M4(u_l/a, 1) + AusmTools::M4(u_r/a, -1);
  real Mp   = 0.5*(M + fabs(M));
  real Mm   = 0.5*(M - fabs(M));
  real p    = AusmTools::P5(u_l/a, 1)*p_l + AusmTools::P5(u_r/a, -1)*p_r;
  countFlops(14);

  flux.var[0] += a*A*(r_l*Mp + r_r*Mm);
  flux.var[1] += 0.5*flux.var[0]*(u_l + u_r) + A*p - 0.5*fabs(flux.var[0])*(u_r - u_l);
  flux.var[2] += 0.5*flux.var[0]*(v_l + v_r)       - 0.5*fabs(flux.var[0])*(v_r - v_l);
  flux.var[3] += 0.5*flux.var[0]*(w_l + w_r)       - 0.5*fabs(flux.var[0])*(w_r - w_l);
  flux.var[4] += 0.5*flux.var[0]*(H_l + H_r)       - 0.5*fabs(flux.var[0])*(H_r - H_l);
  countFlops(36);
}

template <class TReconstruction>
void AusmPlus<TReconstruction>::y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
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

  real a    = 0.5*(a_l + a_r);
  real M    = AusmTools::M4(v_l/a, 1) + AusmTools::M4(v_r/a, -1);
  real Mp   = 0.5*(M + fabs(M));
  real Mm   = 0.5*(M - fabs(M));
  real p    = AusmTools::P5(v_l/a, 1)*p_l + AusmTools::P5(v_r/a, -1)*p_r;
  countFlops(14);

  flux.var[0] += a*A*(r_l*Mp + r_r*Mm);
  flux.var[1] += 0.5*flux.var[0]*(u_l + u_r)       - 0.5*fabs(flux.var[0])*(u_r - u_l);
  flux.var[2] += 0.5*flux.var[0]*(v_l + v_r) + A*p - 0.5*fabs(flux.var[0])*(v_r - v_l);
  flux.var[3] += 0.5*flux.var[0]*(w_l + w_r)       - 0.5*fabs(flux.var[0])*(w_r - w_l);
  flux.var[4] += 0.5*flux.var[0]*(H_l + H_r)       - 0.5*fabs(flux.var[0])*(H_r - H_l);
  countFlops(36);
}

template <class TReconstruction>
void AusmPlus<TReconstruction>::z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
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

  real a    = 0.5*(a_l + a_r);
  real M    = AusmTools::M4(w_l/a, 1) + AusmTools::M4(w_r/a, -1);
  real Mp   = 0.5*(M + fabs(M));
  real Mm   = 0.5*(M - fabs(M));
  real p    = AusmTools::P5(w_l/a, 1)*p_l + AusmTools::P5(w_r/a, -1)*p_r;
  countFlops(14);

  flux.var[0] += a*A*(r_l*Mp + r_r*Mm);
  flux.var[1] += 0.5*flux.var[0]*(u_l + u_r)       - 0.5*fabs(flux.var[0])*(u_r - u_l);
  flux.var[2] += 0.5*flux.var[0]*(v_l + v_r)       - 0.5*fabs(flux.var[0])*(v_r - v_l);
  flux.var[3] += 0.5*flux.var[0]*(w_l + w_r) + A*p - 0.5*fabs(flux.var[0])*(w_r - w_l);
  flux.var[4] += 0.5*flux.var[0]*(H_l + H_r)       - 0.5*fabs(flux.var[0])*(H_r - H_l);
  countFlops(36);
}

#endif // AUSMPLUS_H
