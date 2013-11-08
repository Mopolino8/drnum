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
#ifndef COMPRESSIBLESMAGORINSKYFLUX_H
#define COMPRESSIBLESMAGORINSKYFLUX_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <unsigned int DIM, unsigned int CS, typename TGas>
class CompressibleSmagorinskyFlux : public CompressibleFlux<DIM, TGas>
{

protected: // methods

  /*
  template <typename PATCH>
  static CUDA_DH real velNormGrad(PATCH *patch, size_t i_field, size_t idx)
  {
    dim_t<DIM> dim;
    patch->getVar(dim, i_field, idx, var);
    real scal = var[1]*sf.wnx + var[2]*sf.wny + var[3]*sf.wnz;
    var[1] -= scal*sf.wnx;
    var[2] -= scal*sf.wny;
    var[3] -= scal*sf.wnz;


    size_t i, j, k;
    patch->ijk(idx, i, j, k);

    size_t idx1 = patch->index(i-1, j, k);
    size_t idx2 = patch->index(i+1, j, k);
    size_t idx3 = patch->index(i, j-1, k);
    size_t idx4 = patch->index(i, j+1, k);
    size_t idx5 = patch->index(i, j, k-1);
    size_t idx6 = patch->index(i, j, k+1);

    real wgt1 = real(!patch->isSplitFace(idx, idx1));
    real wgt2 = real(!patch->isSplitFace(idx, idx2));
    real wgt3 = real(!patch->isSplitFace(idx, idx3));
    real wgt4 = real(!patch->isSplitFace(idx, idx4));
    real wgt5 = real(!patch->isSplitFace(idx, idx5));
    real wgt6 = real(!patch->isSplitFace(idx, idx6));

    wgt1 = -wgt1*sf.wnx;
    wgt2 =  wgt2*sf.wnx;
    wgt3 = -wgt3*sf.wny;
    wgt4 =  wgt4*sf.wny;
    wgt5 = -wgt5*sf.wnz;
    wgt6 =  wgt6*sf.wnz;

    if (wgt1 < 0.1) wgt1 = 0;
    if (wgt2 < 0.1) wgt2 = 0;
    if (wgt3 < 0.1) wgt3 = 0;
    if (wgt4 < 0.1) wgt4 = 0;
    if (wgt5 < 0.1) wgt5 = 0;
    if (wgt6 < 0.1) wgt6 = 0;


    real iwgt = wgt1 + wgt2 + wgt3 + wgt4 + wgt5 + wgt6;
    if (iwgt <= 0) return 0;
    iwgt = 1.0/iwgt;

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
      var[i_var] *= iwgt;
    }
  }
  */


public: // methods

  template <typename PATCH>
  CUDA_DH void xField(PATCH *patch,
                      size_t i, size_t j, size_t k,
                      real x, real y, real z,
                      real A, real* flux)
  {
    dim_t<DIM> dim;

    real var_l[DIM], var_r[DIM], var[DIM];
    patch->getVar(dim, 0, i-1, j, k, var_l);
    patch->getVar(dim, 0, i,   j, k, var_r);
    for (int i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = (real)0.5*(var_l[i_var] + var_r[i_var]);
    }
    countFlops(10);

    // Compute some useful coefficients
    real gamma  = TGas::gamma(var);
    real gamma1 = gamma - (real)1.0;
    real ir     = (real)1.0/var[0];
    countFlops(2);

    real gradx[DIM];
    real grady[DIM], grady_l[DIM], grady_r[DIM];
    real gradz[DIM], gradz_l[DIM], gradz_r[DIM];
    patch->getYGrad(dim, 0, i-1, j, k, grady_l);
    patch->getZGrad(dim, 0, i-1, j, k, gradz_l);
    patch->getYGrad(dim, 0, i, j, k, grady_r);
    patch->getZGrad(dim, 0, i, j, k, gradz_r);
    real D = (real)1.0/patch->dx();
    for (int i_var = 0; i_var < DIM; ++i_var) {
      gradx[i_var] = (var_r[i_var] - var_l[i_var])*D;
      grady[i_var] = (real)0.5*(grady_l[i_var] + grady_r[i_var]);
      gradz[i_var] = (real)0.5*(gradz_l[i_var] + gradz_r[i_var]);
    }
    countFlops(31);

    // Compute primitive gradients d(u)/d(.), d(v)/d(.), d(T)/d(.)
    real du_dx = (gradx[1] - gradx[0]*var[1]*ir)*ir;
    real dv_dx = (gradx[2] - gradx[0]*var[2]*ir)*ir;
    real dw_dx = (gradx[3] - gradx[0]*var[3]*ir)*ir;
    real du_dy = (grady[1] - grady[0]*var[1]*ir)*ir;
    real dv_dy = (grady[2] - grady[0]*var[2]*ir)*ir;
    real dw_dy = (grady[3] - grady[0]*var[3]*ir)*ir;
    real du_dz = (gradz[1] - gradz[0]*var[1]*ir)*ir;
    real dv_dz = (gradz[2] - gradz[0]*var[2]*ir)*ir;
    real dw_dz = (gradz[3] - gradz[0]*var[3]*ir)*ir;
    real dT_dx = gamma1*ir*ir*(  gradx[4]*var[0] - gradx[0]*var[4] - var[1]*gradx[1]
                               - var[2]*gradx[2] - var[3]*gradx[3]
                               + gradx[0]*ir*(var[1]*var[1] + var[2]*var[2] + var[3]*var[3]));
    countFlops(56);

    real cs  = 1e-4*real(CS);
    real S   = sqr(du_dx) + sqr(du_dy + dv_dx) + sqr(du_dz + dw_dx) + sqr(dv_dy) + sqr(dv_dz + dw_dy) + sqr(dw_dz);
    real mut = var[0]*sqr(cs*patch->dx())*sqrt(S);

    real F2,F3,F4;
    real mu = TGas::mu(var);
    real Pr = TGas::Pr(var);
    real cp = TGas::cp(var);
    real R  = TGas::R(var);
    real k_heat = mu*cp/Pr + mut*cp/0.85;
    F2 = FR23*(mu + mut);
    F3 = mu + mut;
    F4 = k_heat/R;
    countFlops(28);
    countSqrts(1);

    // tensions
    real tau_xx = F2 * ((real)2.0*du_dx - dv_dy - dw_dz);
    real tau_xy = F3 * ( du_dy + dv_dx);
    real tau_xz = F3 * ( du_dz + dw_dx);
    countFlops(9);

    // heat flux
    real q_dot_x = F4 * dT_dx;
    countFlops(1);

    // The viscous fluxes for the cons. quantities
    flux[1] -= tau_xx*A;
    flux[2] -= tau_xy*A;
    flux[3] -= tau_xz*A;
    flux[4] -= ((tau_xx*var[1] + tau_xy*var[2] + tau_xz*var[3])*ir + q_dot_x)*A;
    countFlops(15);
  }

  template <typename PATCH>
  CUDA_DH void yField(PATCH *patch,
                      size_t i, size_t j, size_t k,
                      real x, real y, real z,
                      real A, real* flux)
  {
    dim_t<DIM> dim;

    real var_l[DIM], var_r[DIM], var[DIM];
    patch->getVar(dim, 0, i, j-1, k, var_l);
    patch->getVar(dim, 0, i, j,   k, var_r);
    for (int i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = (real)0.5*(var_l[i_var] + var_r[i_var]);
    }
    countFlops(10);

    // Compute some useful coefficients
    real gamma  = TGas::gamma(var);
    real gamma1 = gamma - (real)1.0;
    real ir     = (real)1.0/var[0];
    countFlops(2);

    real gradx[DIM], gradx_l[DIM], gradx_r[DIM];
    real grady[DIM];
    real gradz[DIM], gradz_l[DIM], gradz_r[DIM];
    patch->getXGrad(dim, 0, i, j-1, k, gradx_l);
    patch->getZGrad(dim, 0, i, j-1, k, gradz_l);
    patch->getXGrad(dim, 0, i, j, k, gradx_r);
    patch->getZGrad(dim, 0, i, j, k, gradz_r);
    real D = (real)1.0/patch->dy();
    for (int i_var = 0; i_var < DIM; ++i_var) {
      gradx[i_var] = (real)0.5*(gradx_l[i_var] + gradx_r[i_var]);
      grady[i_var] = (var_r[i_var] - var_l[i_var])*D;
      gradz[i_var] = (real)0.5*(gradz_l[i_var] + gradz_r[i_var]);
    }
    countFlops(31);

    // Compute primitive gradients d(u)/d(.), d(v)/d(.), d(T)/d(.)
    real du_dx = (gradx[1] - gradx[0]*var[1]*ir)*ir;
    real dv_dx = (gradx[2] - gradx[0]*var[2]*ir)*ir;
    real dw_dx = (gradx[3] - gradx[0]*var[3]*ir)*ir;
    real du_dy = (grady[1] - grady[0]*var[1]*ir)*ir;
    real dv_dy = (grady[2] - grady[0]*var[2]*ir)*ir;
    real dw_dy = (grady[3] - grady[0]*var[3]*ir)*ir;
    real du_dz = (gradz[1] - gradz[0]*var[1]*ir)*ir;
    real dv_dz = (gradz[2] - gradz[0]*var[2]*ir)*ir;
    real dw_dz = (gradz[3] - gradz[0]*var[3]*ir)*ir;
    real dT_dy = gamma1*ir*ir*(  grady[4]*var[0] - grady[0]*var[4] - var[1]*grady[1]
                                 - var[2]*grady[2] - var[3]*grady[3]
                                 + grady[0]*ir*(var[1]*var[1] + var[2]*var[2] + var[3]*var[3]));
    countFlops(56);

    real cs  = 1e-4*real(CS);
    real S   = sqr(du_dx) + sqr(du_dy + dv_dx) + sqr(du_dz + dw_dx) + sqr(dv_dy) + sqr(dv_dz + dw_dy) + sqr(dw_dz);
    real mut = var[0]*sqr(cs*patch->dy())*sqrt(S);

    real F2,F3,F4;
    real mu = TGas::mu(var);
    real Pr = TGas::Pr(var);
    real cp = TGas::cp(var);
    real R  = TGas::R(var);
    real k_heat = mu*cp/Pr + mut*cp/0.85;
    F2 = FR23*(mu + mut);
    F3 = mu + mut;
    F4 = k_heat/R;
    countFlops(28);
    countSqrts(1);

    // tensions
    real tau_yy = F2 * (-du_dx + (real)2.0*dv_dy - dw_dz);
    real tau_xy = F3 * ( du_dy + dv_dx);
    real tau_yz = F3 * ( dv_dz + dw_dy);
    countFlops(9);

    // heat flux
    real q_dot_y = F4 * dT_dy;
    countFlops(1);

    // The viscous fluxes for the cons. quantities
    flux[1] -= tau_xy*A;
    flux[2] -= tau_yy*A;
    flux[3] -= tau_yz*A;
    flux[4] -= ((tau_xy*var[1] + tau_yy*var[2] + tau_yz*var[3])*ir + q_dot_y)*A;
    countFlops(15);
  }

  template <typename PATCH>
  CUDA_DH void zField(PATCH *patch,
                      size_t i, size_t j, size_t k,
                      real x, real y, real z,
                      real A, real* flux)
  {
    dim_t<DIM> dim;

    real var_l[DIM], var_r[DIM], var[DIM];
    patch->getVar(dim, 0, i, j, k-1, var_l);
    patch->getVar(dim, 0, i, j, k,   var_r);
    for (int i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = (real)0.5*(var_l[i_var] + var_r[i_var]);
    }
    countFlops(10);

    // Compute some useful coefficients
    real gamma  = TGas::gamma(var);
    real gamma1 = gamma - (real)1.0;
    real ir     = (real)1.0/var[0];
    countFlops(2);

    real gradx[DIM], gradx_l[DIM], gradx_r[DIM];
    real grady[DIM], grady_l[DIM], grady_r[DIM];
    real gradz[DIM];
    patch->getXGrad(dim, 0, i, j, k-1, gradx_l);
    patch->getYGrad(dim, 0, i, j, k-1, grady_l);
    patch->getXGrad(dim, 0, i, j, k, gradx_r);
    patch->getYGrad(dim, 0, i, j, k, grady_r);
    real D = (real)1.0/patch->dz();
    for (int i_var = 0; i_var < DIM; ++i_var) {
      gradx[i_var] = (real)0.5*(gradx_l[i_var] + gradx_r[i_var]);
      grady[i_var] = (real)0.5*(grady_l[i_var] + grady_r[i_var]);
      gradz[i_var] = (var_r[i_var] - var_l[i_var])*D;
    }
    countFlops(31);

    // Compute primitive gradients d(u)/d(.), d(v)/d(.), d(T)/d(.)
    real du_dx = (gradx[1] - gradx[0]*var[1]*ir)*ir;
    real dv_dx = (gradx[2] - gradx[0]*var[2]*ir)*ir;
    real dw_dx = (gradx[3] - gradx[0]*var[3]*ir)*ir;
    real du_dy = (grady[1] - grady[0]*var[1]*ir)*ir;
    real dv_dy = (grady[2] - grady[0]*var[2]*ir)*ir;
    real dw_dy = (grady[3] - grady[0]*var[3]*ir)*ir;
    real du_dz = (gradz[1] - gradz[0]*var[1]*ir)*ir;
    real dv_dz = (gradz[2] - gradz[0]*var[2]*ir)*ir;
    real dw_dz = (gradz[3] - gradz[0]*var[3]*ir)*ir;
    real dT_dz = gamma1*ir*ir*(  gradz[4]*var[0] - gradz[0]*var[4] - var[1]*gradz[1]
                                 - var[2]*gradz[2] - var[3]*gradz[3]
                                 + gradz[0]*ir*(var[1]*var[1] + var[2]*var[2] + var[3]*var[3]));
    countFlops(56);

    real cs  = 1e-4*real(CS);
    real S   = sqr(du_dz) + sqr(du_dy + dv_dx) + sqr(du_dz + dw_dx) + sqr(dv_dy) + sqr(dv_dz + dw_dy) + sqr(dw_dz);
    real mut = var[0]*sqr(cs*patch->dx())*sqrt(S);

    real F2,F3,F4;
    real mu = TGas::mu(var);
    real Pr = TGas::Pr(var);
    real cp = TGas::cp(var);
    real R  = TGas::R(var);
    real k_heat = mu*cp/Pr + mut*cp/0.85;
    F2 = FR23*(mu + mut);
    F3 = mu + mut;
    F4 = k_heat/R;
    countFlops(28);
    countSqrts(1);

    // tensions
    real tau_zz = F2 * (-du_dx - dv_dy + (real)2.0*dw_dz);
    real tau_xz = F3 * ( du_dz + dw_dx);
    real tau_yz = F3 * ( dv_dz + dw_dy);
    countFlops(9);

    // heat flux
    real q_dot_z = F4 * dT_dz;
    countFlops(1);

    // The viscous fluxes for the cons. quantities
    flux[1] -= tau_xz*A;
    flux[2] -= tau_yz*A;
    flux[3] -= tau_zz*A;
    flux[4] -= ((tau_xz*var[1] + tau_yz*var[2] + tau_zz*var[3])*ir + q_dot_z)*A;
    countFlops(15);
  }

  template <typename PATCH> CUDA_DH void splitFlux(PATCH *patch, splitface_t sf, real* flux)
  {
    real var[DIM];
    dim_t<DIM> dim;
    patch->getVar(dim, 0, sf.idx, var);
    real u = var[1]/var[0];
    real v = var[2]/var[0];
    real w = var[3]/var[0];
    real u_abs = sqrt(u*u + v*v + w*w);
    real scal = u*sf.wnx + v*sf.wny + w*sf.wnz;
    u -= scal*sf.wnx;
    v -= scal*sf.wny;
    w -= scal*sf.wnz;
    real u_tgt = sqrt(u*u + v*v + w*w);

    if (u_tgt/u_abs < 1e-3 || u_abs < 1e-6) {
      return;
    }

    real y_plus = 20;
    real u_tau = 0;
    real mu = TGas::mu(var);
    real y_wall = sf.wnx*patch->dx() + sf.wny*patch->dy() + sf.wnz*patch->dx();
    for (int iter = 0; iter < 10; ++iter) {
      y_plus = max(20.0, u_tau*y_wall*var[0]/mu);
      u_tau  = u_tgt/(log(y_plus)/0.41 + 5.0);
    }
    real tau_wall = var[0]*u_tau*u_tau*sqrt(sf.nx*sf.nx + sf.ny*sf.ny + sf.nz*sf.nz);
    flux[1] -= u/u_tgt * tau_wall;
    flux[2] -= v/u_tgt * tau_wall;
    flux[3] -= w/u_tgt * tau_wall;
  }


};

#endif // COMPRESSIBLESMAGORINSKYFLUX_H
