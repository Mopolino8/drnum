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
#ifndef COMPRESSIBLEVISCFLUX_H
#define COMPRESSIBLEVISCFLUX_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <typename TGas>
class CompressibleViscFlux : public CompressibleFlux
{

public: // methods

  template <typename PATCH>
  CUDA_DH void xField(PATCH *patch,
                      size_t i, size_t j, size_t k,
                      real x, real y, real z,
                      real A, real* flux)
  {
    real var_l[5], var_r[5], var[5];
    patch->getVar(0, i-1, j, k, var_l);
    patch->getVar(0, i,   j, k, var_r);
    for (int i_var = 0; i_var < 5; ++i_var) {
      var[i_var] = 0.5*(var_l[i_var] + var_r[i_var]);
    }
    countFlops(10);

    // Compute some useful coefficients
    real gamma  = TGas::gamma(var);
    real gamma1 = gamma - 1.0;
    real ir     = 1.0/var[0];
    countFlops(2);

    real F2,F3,F4;
    real mu = TGas::mu(var);
    real Pr = TGas::Pr(var);
    real cp = TGas::cp(var);
    real R  = TGas::R(var);
    real k_heat = mu*cp/Pr;
    F2 = FR23*mu;
    F3 = mu;
    F4 = k_heat/R;
    countFlops(4);

    real gradx[5];
    real grady[5], grady_l[5], grady_r[5];
    real gradz[5], gradz_l[5], gradz_r[5];
    patch->getYGrad(0, i-1, j, k, grady_l);
    patch->getZGrad(0, i-1, j, k, gradz_l);
    patch->getYGrad(0, i, j, k, grady_r);
    patch->getZGrad(0, i, j, k, gradz_r);
    real D = 1.0/patch->dx();
    for (int i_var = 0; i_var < 5; ++i_var) {
      gradx[i_var] = (var_r[i_var] - var_l[i_var])*D;
      grady[i_var] = 0.5*(grady_l[i_var] + grady_r[i_var]);
      gradz[i_var] = 0.5*(gradz_l[i_var] + gradz_r[i_var]);
    }
    countFlops(31);

    // Compute primitive gradients d(u)/d(.), d(v)/d(.), d(T)/d(.)
    real du_dx = (gradx[1] - gradx[0]*var[1]*ir)*ir;
    real dv_dx = (gradx[2] - gradx[0]*var[2]*ir)*ir;
    real dw_dx = (gradx[3] - gradx[0]*var[3]*ir)*ir;
    real dT_dx = gamma1*ir*ir*(  gradx[4]*var[0] - gradx[0]*var[4] - var[1]*gradx[1]
                               - var[2]*gradx[2] - var[3]*gradx[3]
                               + gradx[0]*ir*(var[1]*var[1] + var[2]*var[2] + var[3]*var[3]));
    real du_dy = (grady[1] - grady[0]*var[1]*ir)*ir;
    real dv_dy = (grady[2] - grady[0]*var[2]*ir)*ir;
    real du_dz = (gradz[1] - gradz[0]*var[1]*ir)*ir;
    real dw_dz = (gradz[3] - gradz[0]*var[3]*ir)*ir;
    countFlops(48);

    // tensions
    real tau_xx = F2 * (2.0*du_dx - dv_dy - dw_dz);
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
    real var_l[5], var_r[5], var[5];
    patch->getVar(0, i, j-1, k, var_l);
    patch->getVar(0, i, j,   k, var_r);
    for (int i_var = 0; i_var < 5; ++i_var) {
      var[i_var] = 0.5*(var_l[i_var] + var_r[i_var]);
    }
    countFlops(10);

    // Compute some useful coefficients
    real gamma  = TGas::gamma(var);
    real gamma1 = gamma - 1.0;
    real ir     = 1.0/var[0];
    countFlops(2);

    real F2,F3,F4;
    real mu = TGas::mu(var);
    real Pr = TGas::Pr(var);
    real cp = TGas::cp(var);
    real R  = TGas::R(var);
    real k_heat = mu*cp/Pr;
    F2 = FR23*mu;
    F3 = mu;
    F4 = k_heat/R;
    countFlops(4);

    real gradx[5], gradx_l[5], gradx_r[5];
    real grady[5];
    real gradz[5], gradz_l[5], gradz_r[5];
    patch->getXGrad(0, i, j-1, k, gradx_l);
    patch->getZGrad(0, i, j-1, k, gradz_l);
    patch->getXGrad(0, i, j, k, gradx_r);
    patch->getZGrad(0, i, j, k, gradz_r);
    real D = 1.0/patch->dy();
    for (int i_var = 0; i_var < 5; ++i_var) {
      gradx[i_var] = 0.5*(gradx_l[i_var] + gradx_r[i_var]);
      grady[i_var] = (var_r[i_var] - var_l[i_var])*D;
      gradz[i_var] = 0.5*(gradz_l[i_var] + gradz_r[i_var]);
    }
    countFlops(31);

    // Compute primitive gradients d(u)/d(.), d(v)/d(.), d(T)/d(.)
    real du_dx = (gradx[1] - gradx[0]*var[1]*ir)*ir;
    real dv_dx = (gradx[2] - gradx[0]*var[2]*ir)*ir;
    real du_dy = (grady[1] - grady[0]*var[1]*ir)*ir;
    real dv_dy = (grady[2] - grady[0]*var[2]*ir)*ir;
    real dw_dy = (grady[3] - grady[0]*var[3]*ir)*ir;
    real dT_dy = gamma1*ir*ir*(  grady[4]*var[0] - grady[0]*var[4] - var[1]*grady[1]
                                 - var[2]*grady[2] - var[3]*grady[3]
                                 + grady[0]*ir*(var[1]*var[1] + var[2]*var[2] + var[3]*var[3]));
    real dv_dz = (gradz[2] - gradz[0]*var[2]*ir)*ir;
    real dw_dz = (gradz[3] - gradz[0]*var[3]*ir)*ir;
    countFlops(48);

    // tensions
    real tau_yy = F2 * (-du_dx + 2.0*dv_dy - dw_dz);
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
    real var_l[5], var_r[5], var[5];
    patch->getVar(0, i, j, k-1, var_l);
    patch->getVar(0, i, j, k,   var_r);
    for (int i_var = 0; i_var < 5; ++i_var) {
      var[i_var] = 0.5*(var_l[i_var] + var_r[i_var]);
    }
    countFlops(10);

    // Compute some useful coefficients
    real gamma  = TGas::gamma(var);
    real gamma1 = gamma - 1.0;
    real ir     = 1.0/var[0];
    countFlops(2);

    real F2,F3,F4;
    real mu = TGas::mu(var);
    real Pr = TGas::Pr(var);
    real cp = TGas::cp(var);
    real R  = TGas::R(var);
    real k_heat = mu*cp/Pr;
    F2 = FR23*mu;
    F3 = mu;
    F4 = k_heat/R;
    countFlops(4);

    real gradx[5], gradx_l[5], gradx_r[5];
    real grady[5], grady_l[5], grady_r[5];
    real gradz[5];
    patch->getXGrad(0, i, j, k-1, gradx_l);
    patch->getYGrad(0, i, j, k-1, grady_l);
    patch->getXGrad(0, i, j, k, gradx_r);
    patch->getYGrad(0, i, j, k, grady_r);
    real D = 1.0/patch->dz();
    for (int i_var = 0; i_var < 5; ++i_var) {
      gradx[i_var] = 0.5*(gradx_l[i_var] + gradx_r[i_var]);
      grady[i_var] = 0.5*(grady_l[i_var] + grady_r[i_var]);
      gradz[i_var] = (var_r[i_var] - var_l[i_var])*D;
    }
    countFlops(31);

    // Compute primitive gradients d(u)/d(.), d(v)/d(.), d(T)/d(.)
    real du_dx = (gradx[1] - gradx[0]*var[1]*ir)*ir;
    real dw_dx = (gradx[3] - gradx[0]*var[3]*ir)*ir;
    real dv_dy = (grady[2] - grady[0]*var[2]*ir)*ir;
    real dw_dy = (grady[3] - grady[0]*var[3]*ir)*ir;
    real du_dz = (gradz[1] - gradz[0]*var[1]*ir)*ir;
    real dv_dz = (gradz[2] - gradz[0]*var[2]*ir)*ir;
    real dw_dz = (gradz[3] - gradz[0]*var[3]*ir)*ir;
    real dT_dz = gamma1*ir*ir*(  gradz[4]*var[0] - gradz[0]*var[4] - var[1]*gradz[1]
                                 - var[2]*gradz[2] - var[3]*gradz[3]
                                 + gradz[0]*ir*(var[1]*var[1] + var[2]*var[2] + var[3]*var[3]));
    countFlops(48);

    // tensions
    real tau_zz = F2 * (-du_dx - dv_dy + 2.0*dw_dz);
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

};




#endif // COMPRESSIBLEVISCFLUX_H
