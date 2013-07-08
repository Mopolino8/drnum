// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef COMPRESSIBLEFLUX_H
#define COMPRESSIBLEFLUX_H

#include "blockcfd.h"

#define COMPR_VARS \
  real r  = var[0]; \
  real ir = CHECKED_REAL(real(1.0)/r); \
  real ru = var[1]; \
  real rv = var[2]; \
  real rw = var[3]; \
  real u  = ru*ir; \
  real v  = rv*ir; \
  real w  = rw*ir; \
  real rE = var[4]; \
  real T  = CHECKED_REAL((rE*ir - real(0.5)*(u*u + v*v + w*w))/TGas::cv(var)); \
  real p  = r*TGas::R(var)*T; \
  countFlops(15);

class CompressibleFlux
{

protected: // methods

  static CUDA_DH real M1(real M,real s)
  {
    countFlops(3);
    return CHECKED_REAL(0.5*(M + s*fabs(M)));
  }

  static CUDA_DH real M2(real M,real s)
  {
    countFlops(4);
    return CHECKED_REAL(0.25*s*sqr(M + s));
  }

  static CUDA_DH real M4(real M,real s)
  {
    if (fabs(M) >= 1) {
      return M1(M, s);
    }
    countFlops(4);
    return CHECKED_REAL(M2(M, s)*(1 - 2*s*M2(M, -s)));
  }

  static CUDA_DH real P5(real M,real s)
  {
    if (fabs(M) >= 1) {
      countFlops(1);
      return CHECKED_REAL(M1(M,s)/M);
    }
    countFlops(6);
    return M2(M, s)*((2*s - M) - 3*s*M*M2(M, -s));
  }

};

#define COMPRESSIBLE_LEFT_VARS \
  REGREAL r_l  = var_l[0]; \
  REGREAL ir_l = CHECKED_REAL(1.0/r_l); \
  REGREAL ru_l = var_l[1]; \
  REGREAL rv_l = var_l[2]; \
  REGREAL rw_l = var_l[3]; \
  REGREAL u_l  = ru_l*ir_l; \
  REGREAL v_l  = rv_l*ir_l; \
  REGREAL w_l  = rw_l*ir_l; \
  REGREAL rE_l = var_l[4]; \
  REGREAL T_l  = CHECKED_REAL((rE_l*ir_l - real(0.5)*(u_l*u_l + v_l*v_l + w_l*w_l))/TGas::cv(var_l)); \
  REGREAL p_l  = r_l*TGas::R(var_l)*T_l; \
  REGREAL a_l  = CHECKED_REAL(sqrt(TGas::gamma(var_l)*TGas::R(var_l)*T_l)); \
  REGREAL H_l  = (rE_l + p_l)/r_l; \
  countFlops(19); \
  countSqrts(1);

#define COMPRESSIBLE_RIGHT_VARS \
  REGREAL r_r  = var_r[0]; \
  REGREAL ir_r = CHECKED_REAL(1.0/r_r); \
  REGREAL ru_r = var_r[1]; \
  REGREAL rv_r = var_r[2]; \
  REGREAL rw_r = var_r[3]; \
  REGREAL u_r  = ru_r*ir_r; \
  REGREAL v_r  = rv_r*ir_r; \
  REGREAL w_r  = rw_r*ir_r; \
  REGREAL rE_r = var_r[4]; \
  REGREAL T_r  = CHECKED_REAL((rE_r*ir_r - real(0.5)*(u_r*u_r + v_r*v_r + w_r*w_r))/TGas::cv(var_r)); \
  REGREAL p_r  = r_r*TGas::R(var_r)*T_r; \
  REGREAL a_r  = CHECKED_REAL(sqrt(TGas::gamma(var_r)*TGas::R(var_r)*T_r)); \
  REGREAL H_r  = (rE_r + p_r)/r_r; \
  countFlops(19); \
  countSqrts(1);

#define COMPRESSIBLE_LEFT_PROJX \
  real var_l[5]; \
  m_Reconstruction.project(patch, var_l, 0, 5, i-1, j, k, i, j, k, x - patch->dx(), y, z, x, y, z); \
  COMPRESSIBLE_LEFT_VARS

#define COMPRESSIBLE_RIGHT_PROJX \
  real var_r[5]; \
  m_Reconstruction.project(patch, var_r, 0, 5, i, j, k, i-1, j, k, x, y, z, x - patch->dx(), y, z); \
  COMPRESSIBLE_RIGHT_VARS

#define COMPRESSIBLE_LEFT_PROJY \
  real var_l[5]; \
  m_Reconstruction.project(patch, var_l, 0, 5, i, j-1, k, i, j, k, x, y - patch->dy(), z, x, y, z); \
  COMPRESSIBLE_LEFT_VARS

#define COMPRESSIBLE_RIGHT_PROJY \
  real var_r[5]; \
  m_Reconstruction.project(patch, var_r, 0, 5, i, j, k, i, j-1, k, x, y, z, x, y - patch->dy(), z); \
  COMPRESSIBLE_RIGHT_VARS

#define COMPRESSIBLE_LEFT_PROJZ \
  real var_l[5]; \
  m_Reconstruction.project(patch, var_l, 0, 5, i, j, k-1, i, j, k, x, y, z - patch->dz(), x, y, z); \
  COMPRESSIBLE_LEFT_VARS

#define COMPRESSIBLE_RIGHT_PROJZ \
  real var_r[5]; \
  m_Reconstruction.project(patch, var_r, 0, 5, i, j, k, i, j, k-1, x, y, z, x, y, z - patch->dz()); \
  COMPRESSIBLE_RIGHT_VARS


#endif // COMPRESSIBLEFLUX_H
