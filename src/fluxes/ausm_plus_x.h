/**
 * @file
 * @brief Code-Block: Compute the Cartesian AUSM+ flux in x direction.
 *
 * The flux will be computed between the two cells: (i, j, k) and (i+1, j, k).
 * This code block assumes that the following local variables are available and contain the correct values:
 */

real r_l   = PROJ_UX(r, i, i+1, j, k);
real ir_l  = 1.0/r_l;
real ru_l  = PROJ_UX(ru, i, i+1, j, k);
real rv_l  = PROJ_UX(rv, i, i+1, j, k);
real rw_l  = PROJ_UX(rw, i, i+1, j, k);
real u_l   = ru_l*ir_l;
real v_l   = rv_l*ir_l;
real w_l   = rw_l*ir_l;
real rE_l  = PROJ_UX(rE, i, i+1, j, k);
real T_l   = (rE_l*ir_l - 0.5*(u_l*u_l + v_l*v_l + w_l*w_l))/gasCv();
real p_l   = r_l*gasR()*T_l;
real a_l   = sqrt(gasGamma()*gasR()*T_l);

real r_r   = PROJ_UX(r, i+1, i, j, k);
real ir_r  = 1.0/r_r;
real ru_r  = PROJ_UX(ru, i+1, i, j, k);
real rv_r  = PROJ_UX(rv, i+1, i, j, k);
real rw_r  = PROJ_UX(rw, i+1, i, j, k);
real u_r   = ru_r*ir_r;
real v_r   = rv_r*ir_r;
real w_r   = rw_r*ir_r;
real rE_r  = PROJ_UX(rE, i+1, i, j, k);
real T_r   = (rE_r*ir_r - 0.5*(u_r*u_r + v_r*v_r + w_r*w_r))/gasCv();
real p_r   = r_r*gasR()*T_r;
real a_r   = sqrt(gasGamma()*gasR()*T_r);

real a  = 0.5*(a_l + a_r);
real M  = AusmTools::M4(u_l/a, 1) + AusmTools::M4(u_r/a, -1);
real Mp = 0.5*(M + fabs(M));
real Mm = 0.5*(M - fabs(M));
real p  = AusmTools::P5(u_l/a, 1)*p_l + AusmTools::P5(u_r/a, -1)*p_r;
real A  = dy()*dz();

real flux_r  = a*A*(r_l*Mp + r_r*Mm);
real flux_ru = 0.5*flux_r*(ru_l + ru_r) + A*p - 0.5*fabs(flux_r)*(ru_r - ru_l);
real flux_rv = 0.5*flux_r*(rv_l + rv_r)       - 0.5*fabs(flux_r)*(rv_r - rv_l);
real flux_rw = 0.5*flux_r*(rw_l + rw_r)       - 0.5*fabs(flux_r)*(rw_r - rw_l);
real flux_rE = 0.5*flux_r*(rE_l + rE_r)       - 0.5*fabs(flux_r)*(rE_r - rE_l);

#ifdef DEBUG
  CHECK_COMPR_FLUX;
#endif

ADD_COMPR_XFLUX;
