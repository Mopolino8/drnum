#ifndef AUSMPLUS_H
#define AUSMPLUS_H

#include "fluxes/ausmbase.h"
#include "cartesianpatch.h"

template <class TReconstruction>
class AusmPlus : public AusmBase
{

public: // methods

  void x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, flux_t &flux);

};


template <class TReconstruction>
void AusmPlus<TReconstruction>::x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, flux_t &flux)
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

  real a    = 0.5*(a_l + a_r);
  real M    = AusmTools::M4(u_l/a, 1) + AusmTools::M4(u_r/a, -1);
  real Mp   = 0.5*(M + fabs(M));
  real Mm   = 0.5*(M - fabs(M));
  real p    = AusmTools::P5(u_l/a, 1)*p_l + AusmTools::P5(u_r/a, -1)*p_r;

  flux.rho  = a*A*(r_l*Mp + r_r*Mm);
  flux.rhou = 0.5*flux.rho*(ru_l + ru_r) + A*p - 0.5*fabs(flux.rho)*(ru_r - ru_l);
  flux.rhov = 0.5*flux.rho*(rv_l + rv_r)       - 0.5*fabs(flux.rho)*(rv_r - rv_l);
  flux.rhow = 0.5*flux.rho*(rw_l + rw_r)       - 0.5*fabs(flux.rho)*(rw_r - rw_l);
  flux.rhoE = 0.5*flux.rho*(rE_l + rE_r)       - 0.5*fabs(flux.rho)*(rE_r - rE_l);
}

#endif // AUSMPLUS_H
