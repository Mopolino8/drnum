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
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef COMPRESSIBLEFLUX_H
#define COMPRESSIBLEFLUX_H

#include "drnum.h"

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

template <unsigned int DIM, typename TGas>
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

  template <typename PATCH> static CUDA_DH bool averageVar(PATCH *patch, size_t i_field, size_t idx, size_t idx_neigh, real *var)
  {
    dim_t<DIM> dim;
    //patch->getVar(dim, i_field, idx, var);

    size_t i, j, k;
    patch->ijk(idx_neigh, i, j, k);

    size_t idx1 = patch->index(i-1, j, k);
    size_t idx2 = patch->index(i+1, j, k);
    size_t idx3 = patch->index(i, j-1, k);
    size_t idx4 = patch->index(i, j+1, k);
    size_t idx5 = patch->index(i, j, k-1);
    size_t idx6 = patch->index(i, j, k+1);
    real   wgt1 = real(patch->isSplitFace(idx_neigh, idx1));
    real   wgt2 = real(patch->isSplitFace(idx_neigh, idx2));
    real   wgt3 = real(patch->isSplitFace(idx_neigh, idx3));
    real   wgt4 = real(patch->isSplitFace(idx_neigh, idx4));
    real   wgt5 = real(patch->isSplitFace(idx_neigh, idx5));
    real   wgt6 = real(patch->isSplitFace(idx_neigh, idx6));

    real wgt = wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6;

    if (wgt < 0.5) {
      return false;
    }

    real var1[DIM];
    real var2[DIM];
    real var3[DIM];
    real var4[DIM];
    real var5[DIM];
    real var6[DIM];

    patch->getVar(dim, i_field, idx1, var1);
    patch->getVar(dim, i_field, idx2, var2);
    patch->getVar(dim, i_field, idx3, var3);
    patch->getVar(dim, i_field, idx4, var4);
    patch->getVar(dim, i_field, idx5, var5);
    patch->getVar(dim, i_field, idx6, var6);

    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      var[i_var]  = wgt1*var1[i_var];
      var[i_var] += wgt2*var2[i_var];
      var[i_var] += wgt3*var3[i_var];
      var[i_var] += wgt4*var4[i_var];
      var[i_var] += wgt5*var5[i_var];
      var[i_var] += wgt6*var6[i_var];
      var[i_var] /= wgt;
    }
    return true;
  }

  template <typename PATCH> static CUDA_DH bool fieldVar(PATCH *patch, size_t i_field, size_t idx, size_t idx_neigh, real *var)
  {
    dim_t<DIM> dim;

    size_t i1, j1, k1;
    size_t i2, j2, k2;
    patch->ijk(idx, i1, j1, k1);
    patch->ijk(idx_neigh, i2, j2, k2);
    size_t i0 = 2*i1 - i2;
    size_t j0 = 2*j1 - j2;
    size_t k0 = 2*k1 - k2;
    if (!patch->checkRange(i0, j0, k0)) {
      i0 = i1;
      j0 = j1;
      k0 = k1;
    }
    size_t idx0 = patch->index(i0, j0, k0);
    if (patch->isSplitFace(idx0, idx)) {
      i0 = i1;
      j0 = j1;
      k0 = k1;
      idx0 = idx;
    }

    patch->getVar(dim, i_field, idx0, var);
  }


public: // methods

  template <typename PATCH> CUDA_DH void splitFlux(PATCH *patch, splitface_t sf, real* flux)
  {
    /*
    real var_field[DIM];
    real var_cell[DIM];
    real var_face[DIM];
    dim_t<DIM> dim;
    fieldVar(patch, 0, sf.idx, sf.idx_neigh, var_field);
    patch->getVar(dim, 0, sf.idx, var_cell);

    real p_field, u_field, v_field, w_field, T_field;
    real p_cell, u_cell, v_cell, w_cell, T_cell;
    TGas::conservativeToPrimitive(var_field, p_field, T_field, u_field, v_field, w_field);
    TGas::conservativeToPrimitive(var_cell, p_cell, T_cell, u_cell, v_cell, w_cell);

    //real u_abs_field  = sqrt(sqr(u_field) + sqr(v_field) + sqr(w_field));
    //real u_abs_cell   = sqrt(sqr(u_cell) + sqr(v_cell) + sqr(w_cell));
    //real u_norm_field = u_field*sf.wnx + v_field*sf.wny + w_field*sf.wnz;
    //real u_norm_cell  = u_cell*sf.wnx + v_cell*sf.wny + w_cell*sf.wnz;
    //real u_tang_field = sqrt(sqr(u_abs_field) - sqr(u_norm_field));
    //real u_tang_cell  = sqrt(sqr(u_abs_cell) - sqr(u_norm_cell));

    real A = sqrt(sqr(sf.nx) + sqr(sf.ny) + sqr(sf.nz));
    real nx0 = sf.nx/A;
    real ny0 = sf.ny/A;
    real nz0 = sf.nz/A;
    real wgt  = 0.0;

    real Fwgt = 1.0;
    real F_field = var_field[0]*sqr(u_field*nx0 + v_field*ny0 + w_field*nz0) + p_field;
    real F_cell = var_cell[0]*sqr(u_cell*nx0 + v_cell*ny0 + w_cell*nz0) + p_cell;
    real F = Fwgt*F_cell + (1.0 - Fwgt)*F_field;

    //real u_norm_face = FR13*(1.0 - sf.h1)*u_norm_field;
    real p_face = wgt*p_field + (1.0 - wgt)*p_cell;
    real T_face = wgt*T_field + (1.0 - wgt)*T_cell;
    real u_face = wgt*u_field + (1.0 - wgt)*u_cell;
    real v_face = wgt*v_field + (1.0 - wgt)*v_cell;
    real w_face = wgt*w_field + (1.0 - wgt)*w_cell;
    real scal = u_face*sf.wnx + v_face*sf.wny + w_face*sf.wnz;
    u_face -= scal*sf.wnx;
    v_face -= scal*sf.wny;
    w_face -= scal*sf.wnz;
    real un_face = u_face*nx0 + v_face*ny0 + w_face*nz0;
    real RT = TGas::R(var_cell)*T_face;
    real pwgt = fabs(nx0*sf.wnx + ny0*sf.wny + nz0*sf.wnz);
    p_face = pwgt*p_face + (1.0 - pwgt)*F*RT/(sqr(un_face) + RT);
    //u_face += u_norm_face*sf.wnx;
    //v_face += u_norm_face*sf.wny;
    //w_face += u_norm_face*sf.wnz;

    TGas::primitiveToConservative(p_face, T_face, u_face, v_face, w_face, var_face);

    real r_face  = var_face[0];
    real rE_face = var_face[4];

    real h_face = (rE_face + p_face)/r_face;
    real flux_0;

    flux_0   = r_face*u_face*sf.nx + r_face*v_face*sf.ny + r_face*w_face*sf.nz;
    flux[0] += flux_0;
    flux[1] += (r_face*u_face*r_face*u_face + p_face)*sf.nx +                   r_face*u_face*v_face*sf.ny +                   r_face*u_face*w_face*sf.nz;
    flux[2] +=                   r_face*u_face*v_face*sf.nx + (r_face*v_face*r_face*v_face + p_face)*sf.ny +                   r_face*v_face*w_face*sf.nz;
    flux[3] +=                   r_face*u_face*w_face*sf.nx +                   r_face*v_face*w_face*sf.ny + (r_face*w_face*r_face*w_face + p_face)*sf.nz;
    flux[4] += h_face*flux_0;
    */
    real var_field[DIM];
    real var_cell[DIM];
    real var_face[DIM];
    dim_t<DIM> dim;
    fieldVar(patch, 0, sf.idx, sf.idx_neigh, var_field);
    patch->getVar(dim, 0, sf.idx, var_cell);
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
  real var_l[DIM]; \
  m_Reconstruction.project(patch, var_l, 0, i-1, j, k, i, j, k, x - patch->dx(), y, z, x, y, z); \
  COMPRESSIBLE_LEFT_VARS

#define COMPRESSIBLE_RIGHT_PROJX \
  real var_r[DIM]; \
  m_Reconstruction.project(patch, var_r, 0, i, j, k, i-1, j, k, x, y, z, x - patch->dx(), y, z); \
  COMPRESSIBLE_RIGHT_VARS

#define COMPRESSIBLE_LEFT_PROJY \
  real var_l[DIM]; \
  m_Reconstruction.project(patch, var_l, 0, i, j-1, k, i, j, k, x, y - patch->dy(), z, x, y, z); \
  COMPRESSIBLE_LEFT_VARS

#define COMPRESSIBLE_RIGHT_PROJY \
  real var_r[DIM]; \
  m_Reconstruction.project(patch, var_r, 0, i, j, k, i, j-1, k, x, y, z, x, y - patch->dy(), z); \
  COMPRESSIBLE_RIGHT_VARS

#define COMPRESSIBLE_LEFT_PROJZ \
  real var_l[DIM]; \
  m_Reconstruction.project(patch, var_l, 0, i, j, k-1, i, j, k, x, y, z - patch->dz(), x, y, z); \
  COMPRESSIBLE_LEFT_VARS

#define COMPRESSIBLE_RIGHT_PROJZ \
  real var_r[DIM]; \
  m_Reconstruction.project(patch, var_r, 0, i, j, k, i, j, k-1, x, y, z, x, y, z - patch->dz()); \
  COMPRESSIBLE_RIGHT_VARS


#endif // COMPRESSIBLEFLUX_H
