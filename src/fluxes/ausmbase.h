#ifndef AUSMBASE_H
#define AUSMBASE_H

#include "compressibleobject.h"

class AusmBase : public CompressibleObject
{

protected: // methods

  real M1(real M,real s);
  real M2(real M,real s);
  real M4(real M,real s);
  real P5(real M,real s);

};

inline real AusmBase::M1(real M,real s)
{
  return 0.5*(M + s*fabs(M));
}

inline real AusmBase::M2(real M,real s)
{
  return 0.25*s*sqr(M + s);
}

inline real AusmBase::M4(real M,real s)
{
  if (fabs(M) >= 1) {
    return M1(M, s);
  }
  return M2(M, s)*(1 - 2*s*M2(M, -s));
}

inline real AusmBase::P5(real M, real s)
{
  if (fabs(M) >= 1) {
    return M1(M,s)/M;
  }
  return M2(M, s)*((2*s - M) - 3*s*M2(M, -s));
}

#define AUSM_LEFT_VARS \
  REGREAL r_l  = var_l[0]; \
  REGREAL ir_l = CHECKED_REAL(1.0/r_l); \
  REGREAL ru_l = var_l[1]; \
  REGREAL rv_l = var_l[2]; \
  REGREAL rw_l = var_l[3]; \
  REGREAL u_l  = ru_l*ir_l; \
  REGREAL v_l  = rv_l*ir_l; \
  REGREAL w_l  = rw_l*ir_l; \
  REGREAL rE_l = var_l[4]; \
  REGREAL T_l  = CHECKED_REAL((rE_l*ir_l - 0.5*(u_l*u_l + v_l*v_l + w_l*w_l))/gasCv()); \
  REGREAL p_l  = r_l*gasR()*T_l; \
  REGREAL a_l  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_l)); \
  REGREAL H_l  = (rE_l + p_l)/r_l; \
  countFlops(19); \
  countSqrts(1);

#define AUSM_RIGHT_VARS \
  REGREAL r_r  = var_r[0]; \
  REGREAL ir_r = CHECKED_REAL(1.0/r_r); \
  REGREAL ru_r = var_r[1]; \
  REGREAL rv_r = var_r[2]; \
  REGREAL rw_r = var_r[3]; \
  REGREAL u_r  = ru_r*ir_r; \
  REGREAL v_r  = rv_r*ir_r; \
  REGREAL w_r  = rw_r*ir_r; \
  REGREAL rE_r = var_r[4]; \
  REGREAL T_r  = CHECKED_REAL((rE_r*ir_r - 0.5*(u_r*u_r + v_r*v_r + w_r*w_r))/gasCv()); \
  REGREAL p_r  = r_r*gasR()*T_r; \
  REGREAL a_r  = CHECKED_REAL(sqrt(gasGamma()*gasR()*T_r)); \
  REGREAL H_r  = (rE_r + p_r)/r_r; \
  countFlops(19); \
  countSqrts(1);

#define AUSM_LEFT_PROJX \
  real var_l[5]; \
  TReconstruction::project(patch, var_l, 0, 5, i-1, j, k, i, j, k); \
  AUSM_LEFT_VARS

#define AUSM_RIGHT_PROJX \
  real var_r[5]; \
  TReconstruction::project(patch, var_r, 0, 5, i, j, k, i-1, j, k); \
  AUSM_RIGHT_VARS

#define AUSM_LEFT_PROJY \
  real var_l[5]; \
  TReconstruction::project(patch, var_l, 0, 5, i, j-1, k, i, j, k); \
  AUSM_LEFT_VARS

#define AUSM_RIGHT_PROJY \
  real var_r[5]; \
  TReconstruction::project(patch, var_r, 0, 5, i, j, k, i, j-1, k); \
  AUSM_RIGHT_VARS

#define AUSM_LEFT_PROJZ \
  real var_l[5]; \
  TReconstruction::project(patch, var_l, 0, 5, i, j, k-1, i, j, k); \
  AUSM_LEFT_VARS

#define AUSM_RIGHT_PROJZ \
  real var_r[5]; \
  TReconstruction::project(patch, var_r, 0, 5, i, j, k, i, j, k-1); \
  AUSM_RIGHT_VARS


#endif // AUSMBASE_H
