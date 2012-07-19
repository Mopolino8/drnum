#ifndef AUSM_H
#define AUSM_H

#include "ausmbase.h"

template <class TReconstruction>
class Ausm : public AusmBase
{

public: // methods

  void x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, flux_t &flux);

};


template <class TReconstruction>
void Ausm<TReconstruction>::x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, flux_t &flux)
{
  real r_l  = TReconstruction::xp(P, m_Rho, i, j, k);
  real ir_l = 1.0/r_l;
  real ru_l = TReconstruction::xp(P, m_Rhou, i, j, k);
  real rv_l = TReconstruction::xp(P, m_Rhov, i, j, k);
  real rw_l = TReconstruction::xp(P, m_Rhow, i, j, k);
  real u_l  = ru_l*ir_l;
  real v_l  = rv_l*ir_l;
  real w_l  = rw_l*ir_l;
  real rE_l = TReconstruction::xp(P, m_RhoE, i, j, k);
  real T_l  = (rE_l*ir_l - 0.5*(u_l*u_l + v_l*v_l + w_l*w_l))/gasCv();
  real p_l  = r_l*gasR()*T_l;
  real a_l  = sqrt(gasGamma()*gasR()*T_l);

  real r_r  = TReconstruction::xm(P, m_Rho, i, j, k);
  real ir_r = 1.0/r_r;
  real ru_r = TReconstruction::xm(P, m_Rhou, i, j, k);
  real rv_r = TReconstruction::xm(P, m_Rhov, i, j, k);
  real rw_r = TReconstruction::xm(P, m_Rhow, i, j, k);
  real u_r  = ru_r*ir_r;
  real v_r  = rv_r*ir_r;
  real w_r  = rw_r*ir_r;
  real rE_r = TReconstruction::xm(P, m_RhoE, i, j, k);
  real T_r  = (rE_r*ir_r - 0.5*(u_r*u_r + v_r*v_r + w_r*w_r))/gasCv();
  real p_r  = r_r*gasR()*T_r;
  real a_r  = sqrt(gasGamma()*gasR()*T_r);

  real Ma      = 0.5*(u_l/a_l + u_r/a_r);
  real alpha_l = max(0.0, min(1.0, 0.5*(1 + Ma)));
  real alpha_r = 1 - alpha_l;
  real p       = alpha_l*p_l + alpha_r*p_r;
  real conv    = 0.5*(u_l + u_r)*A;

  if (Ma > 0) {
    flux.rho  = conv*r_l;
    flux.rhou = conv*ru_l + A*p;
    flux.rhov = conv*rv_l;
    flux.rhow = conv*rw_l;
    flux.rhoE = conv*(rE_l + p_l);
  } else {
    flux.rho  = conv*r_r;
    flux.rhou = conv*ru_r + A*p;
    flux.rhov = conv*rv_r;
    flux.rhow = conv*rw_r;
    flux.rhoE = conv*(rE_r + p_r);
  }
}

#endif // AUSM_H
