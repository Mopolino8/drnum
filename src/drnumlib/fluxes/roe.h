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
#ifndef ROE_H
#define ROE_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <unsigned int DIM, typename TReconstruction, typename TGas>
class Roe : public CompressibleFlux<DIM, TGas>
{

  TReconstruction m_Reconstruction;

public: // methods

  CUDA_DH real correctLambda(real L)
  {
    const real delta = 0.1;
    L = fabs(L);
    if (L < delta) {
      countFlops(3);
      return 0.5*(sqr(L)/delta + delta);
    }
    return L;
  }

  template <typename PATCH> CUDA_DH void xField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJX;
    COMPRESSIBLE_RIGHT_PROJX;

    real alpha1, alpha2, alpha3, alpha4, alpha5, alpha;
    real a1, a2, a3, a4, a5, b1, b2, b3, b4, b5;
    real ep_r, ep_l;
    real ruv_r, ruv_l, ruw_r, ruw_l, rvw_r, rvw_l;
    real lambda_l, lambda_r;
    real delta_r, delta_ru, delta_rv, delta_rw, delta_rE;
    real hnx;
    real tempo11, usc;

    real flux_r  = 0;
    real flux_ru = 0;
    real flux_rv = 0;
    real flux_rw = 0;
    real flux_rE = 0;

    real nx = A;

    real gam_l = TGas::gamma(var_l);
    real gam_r = TGas::gamma(var_r);
    real gam   = 0.5*(gam_l + gam_r);
    real gam1  = gam - 1.0;

    //
    //.. Compute the ROE Averages
    //
    //.... some useful quantities
    real coef1 = sqrt(r_l);
    real coef2 = sqrt(r_r);
    real somme_coef = coef1 + coef2;
    real isomme_coef = 1.0/somme_coef;
    real h_l = (gam * rE_l*ir_l) - (.5*gam1) * (u_l * u_l + v_l * v_l + w_l * w_l);
    real h_r = (gam * rE_r*ir_r) - (.5*gam1) * (u_r * u_r + v_r * v_r + w_r * w_r);

    //.... averages
    //real rho_ave = coef1 * coef2;
    real u_ave = (coef1 * u_l + coef2 * u_r) * isomme_coef;
    real v_ave = (coef1 * v_l + coef2 * v_r) * isomme_coef;
    real w_ave = (coef1 * w_l + coef2 * w_r) * isomme_coef;
    real h_ave = (coef1 * h_l + coef2 * h_r) * isomme_coef;

    //
    //.. Compute Speed of sound
    real scal     = u_ave*nx;
    real iA       = 1.0/A;
    real u2pv2pw2 = u_ave*u_ave + v_ave*v_ave + w_ave*w_ave;
    real c_speed  = gam1 * (h_ave - 0.5 * u2pv2pw2);
    if(c_speed < 1e-6) {
      c_speed = 1e-6; // avoid division by 0 if critical
    }
    c_speed = sqrt(c_speed);
    real c_speed2 = c_speed*A;

    //
    //.. Compute the eigenvalues of the Jacobian matrix
    real eig_val1 = scal - c_speed2;
    real eig_val2 = scal;
    real eig_val3 = scal + c_speed2;
    //
    //.. Compute the ROE flux
    //.... In this part many tests upon the eigenvalues
    //.... are done to simplify calculations
    //.... Here we use the two formes of the ROE flux :
    //.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
    //.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
    //
    if (eig_val2 <= 0.0) {
      p_r = gam1*(rE_r - 0.5*(ru_r*ru_r + rv_r*rv_r + rw_r*rw_r)*ir_r);
      ep_r = rE_r + p_r;
      ruv_r = ru_r*v_r;
      ruw_r = ru_r*w_r;
      rvw_r = rv_r*w_r;
      flux_r  = ru_r*nx;
      flux_ru = (ru_r*u_r + p_r)*nx;
      flux_rv = ruv_r*nx;
      flux_rw = ruw_r*nx;
      flux_rE = ep_r*u_r*nx;
      //
      //.... A Entropic modification
      //
      p_l = gam1*(rE_l - 0.5*(ru_l*ru_l + rv_l*rv_l + rw_l*rw_l)*ir_l);
      lambda_l = u_l*nx + A*sqrt(gam*p_l*ir_l);
      lambda_r = u_r*nx + A*sqrt(gam*p_r*ir_r);
      if ((lambda_l < 0.) && (lambda_r > 0.)) {
        eig_val3 = lambda_r * (eig_val3 - lambda_l) / (lambda_r - lambda_l);
      };
      //
      if (eig_val3 > 0.0) {
        //.. In this case A+ is obtained by multiplying the last
        //.. colomne of T-1 with the last row of T with eig_val3
        delta_r  = r_r  - r_l;
        delta_ru = ru_r - ru_l;
        delta_rv = rv_r - rv_l;
        delta_rw = rw_r - rw_l;
        delta_rE = rE_r - rE_l;
        //
        scal = scal * iA;
        hnx = nx*iA;
        usc = 1.0/c_speed;
        tempo11 = gam1*usc;
        //.. Last columne of the matrix T-1
        a1 = usc;
        a2 = u_ave*usc + hnx;
        a3 = v_ave*usc;
        a4 = w_ave*usc;
        a5 = 0.5*u2pv2pw2*usc + 2.5*c_speed + scal;
        //.. Last row of the matrix T * eig_val3
        b1 = 0.5*(0.5*tempo11*u2pv2pw2 - scal);
        b2 = 0.5*(hnx - tempo11*u_ave);
        b3 = -0.5*tempo11*v_ave;
        b4 = -0.5*tempo11*w_ave;
        b5 = 0.5*tempo11;
        //
        alpha1 = b1*delta_r;
        alpha2 = b2*delta_ru;
        alpha3 = b3*delta_rv;
        alpha4 = b4*delta_rw;
        alpha5 = b5*delta_rE;
        alpha  = eig_val3*(alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
        //
        flux_r  -= a1*alpha;
        flux_ru -= a2*alpha;
        flux_rv -= a3*alpha;
        flux_rw -= a4*alpha;
        flux_rE -= a5*alpha;
      }
    }
    //
    if(eig_val2 > 0.0) {
      p_l = gam1*(rE_l - 0.5*(ru_l*ru_l + rv_l*rv_l + rw_l*rw_l)*ir_l);
      ep_l = rE_l + p_l;
      ruv_l = ru_l*v_l;
      ruw_l = ru_l*w_l;
      rvw_l = rv_l*w_l;
      flux_r  = ru_l*nx;
      flux_ru = (ru_l*u_l + p_l)*nx;
      flux_rv = ruv_l*nx;
      flux_rw = ruw_l*nx;
      flux_rE = ep_l*u_l*nx;
      //
      // A Entropic modification
      //
      p_r = gam1*(rE_r - 0.5*(ru_r*ru_r + rv_r*rv_r + rw_r*rw_r)*ir_r);
      lambda_l = u_l*nx - A*sqrt(gam*p_l*ir_l);
      lambda_r = u_r*nx - A*sqrt(gam*p_r*ir_r);
      if ((lambda_l < 0.) && (lambda_r > 0.)) {
        eig_val1 = lambda_l * (lambda_r - eig_val1) / (lambda_r - lambda_l);
      };
      //
      if (eig_val1 < 0.0) {
        //.. In this case A+ is obtained by multiplying the first
        //.. columne of T-1 with the first row of T with eig_val1
        delta_r  = r_r  - r_l;
        delta_ru = ru_r - ru_l;
        delta_rv = rv_r - rv_l;
        delta_rw = rw_r - rw_l;
        delta_rE = rE_r - rE_l;
        //
        scal = scal * iA;
        hnx = nx * iA;
        usc = 1.0/c_speed;
        tempo11 = gam1 * usc;
        //.. First colomne of the matrix T-1
        a1 = usc;
        a2 = u_ave * usc - hnx;
        a3 = v_ave * usc;
        a4 = w_ave * usc;
        a5 = 0.5 * u2pv2pw2 * usc + 2.5 * c_speed - scal;
        //.. First row of the matrix T * eig_val1
        b1 = 0.5 * (0.5 * tempo11 * u2pv2pw2 + scal);
        b2 = -0.5 * (hnx + tempo11 * u_ave);
        b3 = -0.5 * (tempo11 * v_ave);
        b4 = -0.5 * (tempo11 * w_ave);
        b5 = 0.5 * tempo11;
        //
        alpha1 = b1 * delta_r;
        alpha2 = b2 * delta_ru;
        alpha3 = b3 * delta_rv;
        alpha4 = b4 * delta_rw;
        alpha5 = b5 * delta_rE;
        alpha  = eig_val1 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
        //
        flux_r  += a1 * alpha;
        flux_ru += a2 * alpha;
        flux_rv += a3 * alpha;
        flux_rw += a4 * alpha;
        flux_rE += a5 * alpha;
      }
    }

    // Add to interface residuals
    flux[0]  += flux_r;
    flux[1] += flux_ru;
    flux[2] += flux_rv;
    flux[3] += flux_rw;
    flux[4] += flux_rE;
  }

  template <typename PATCH> CUDA_DH void yField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJY;
    COMPRESSIBLE_RIGHT_PROJY;

    real alpha1, alpha2, alpha3, alpha4, alpha5, alpha;
    real a1, a2, a3, a4, a5, b1, b2, b3, b4, b5;
    real ep_r, ep_l;
    real ruv_r, ruv_l, ruw_r, ruw_l, rvw_r, rvw_l;
    real lambda_l, lambda_r;
    real delta_r, delta_ru, delta_rv, delta_rw, delta_rE;
    real hny;
    real tempo11, usc;

    real flux_r  = 0;
    real flux_ru = 0;
    real flux_rv = 0;
    real flux_rw = 0;
    real flux_rE = 0;

    real ny = A;

    real gam_l = TGas::gamma(var_l);
    real gam_r = TGas::gamma(var_r);
    real gam   = 0.5*(gam_l + gam_r);
    real gam1  = gam - 1.0;

    //
    //.. Compute the ROE Averages
    //
    //.... some useful quantities
    real coef1 = sqrt(r_l);
    real coef2 = sqrt(r_r);
    real somme_coef = coef1 + coef2;
    real isomme_coef = 1.0/somme_coef;
    real h_l = (gam * rE_l*ir_l) - (.5*gam1) * (u_l * u_l + v_l * v_l + w_l * w_l);
    real h_r = (gam * rE_r*ir_r) - (.5*gam1) * (u_r * u_r + v_r * v_r + w_r * w_r);

    //.... averages
    //real rho_ave = coef1 * coef2;
    real u_ave = (coef1 * u_l + coef2 * u_r) * isomme_coef;
    real v_ave = (coef1 * v_l + coef2 * v_r) * isomme_coef;
    real w_ave = (coef1 * w_l + coef2 * w_r) * isomme_coef;
    real h_ave = (coef1 * h_l + coef2 * h_r) * isomme_coef;

    //
    //.. Compute Speed of sound
    real scal = v_ave * ny;
    real iA = 1.0/A;
    real u2pv2pw2 = u_ave * u_ave + v_ave * v_ave + w_ave * w_ave;
    real c_speed = gam1 * (h_ave - 0.5 * u2pv2pw2);
    if(c_speed < 1e-6) c_speed = 1e-6;    // avoid division by 0 if critical
    c_speed = sqrt(c_speed);
    real c_speed2 = c_speed * A;

    //
    //.. Compute the eigenvalues of the Jacobian matrix
    real eig_val1 = scal - c_speed2;
    real eig_val2 = scal;
    real eig_val3 = scal + c_speed2;
    //
    //.. Compute the ROE flux
    //.... In this part many tests upon the eigenvalues
    //.... are done to simplify calculations
    //.... Here we use the two formes of the ROE flux :
    //.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
    //.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
    //
    if (eig_val2 <= 0.0) {
      p_r = gam1 * (rE_r - 0.5 * (ru_r * ru_r + rv_r * rv_r + rw_r * rw_r) * ir_r);
      ep_r = rE_r + p_r;
      ruv_r = ru_r * v_r;
      ruw_r = ru_r * w_r;
      rvw_r = rv_r * w_r;
      flux_r  = rv_r * ny;
      flux_ru = ruv_r * ny;
      flux_rv = (rv_r * v_r + p_r) * ny;
      flux_rw = rvw_r * ny;
      flux_rE = ep_r * v_r * ny;
      //
      //.... A Entropic modification
      //
      p_l = gam1 * (rE_l - 0.5 * (ru_l * ru_l + rv_l * rv_l + rw_l * rw_l) * ir_l);
      lambda_l = v_l * ny + A * sqrt(gam * p_l * ir_l);
      lambda_r = v_r * ny + A * sqrt(gam * p_r * ir_r);
      if ((lambda_l < 0.) && (lambda_r > 0.)) {
        eig_val3 = lambda_r * (eig_val3 - lambda_l) / (lambda_r - lambda_l);
      };
      //
      if (eig_val3 > 0.0) {
        //.. In this case A+ is obtained by multiplying the last
        //.. colomne of T-1 with the last row of T with eig_val3
        delta_r  = r_r  - r_l;
        delta_ru = ru_r - ru_l;
        delta_rv = rv_r - rv_l;
        delta_rw = rw_r - rw_l;
        delta_rE = rE_r - rE_l;
        //
        scal = scal * iA;
        hny = ny * iA;
        usc = 1.0/c_speed;
        tempo11 = gam1 * usc;
        //.. Last columne of the matrix T-1
        a1 = usc;
        a2 = u_ave * usc;
        a3 = v_ave * usc + hny;
        a4 = w_ave * usc;
        a5 = 0.5 * u2pv2pw2 * usc + 2.5 * c_speed + scal;
        //.. Last row of the matrix T * eig_val3
        b1 = 0.5 * (0.5 * tempo11 * u2pv2pw2 - scal);
        b2 = -0.5 * tempo11 * u_ave;
        b3 = 0.5 * (hny - tempo11 * v_ave);
        b4 = -0.5 * tempo11 * w_ave;
        b5 = 0.5 * tempo11;
        //
        alpha1 = b1 * delta_r;
        alpha2 = b2 * delta_ru;
        alpha3 = b3 * delta_rv;
        alpha4 = b4 * delta_rw;
        alpha5 = b5 * delta_rE;
        alpha  = eig_val3 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
        //
        flux_r  -= a1 * alpha;
        flux_ru -= a2 * alpha;
        flux_rv -= a3 * alpha;
        flux_rw -= a4 * alpha;
        flux_rE -= a5 * alpha;
      };
    };
    //
    if(eig_val2 > 0.0) {
      p_l = gam1 * (rE_l - 0.5 * (ru_l * ru_l + rv_l * rv_l + rw_l * rw_l) * ir_l);
      ep_l = rE_l + p_l;
      ruv_l = ru_l * v_l;
      ruw_l = ru_l * w_l;
      rvw_l = rv_l * w_l;
      flux_r  = rv_l * ny ;
      flux_ru = ruv_l * ny;
      flux_rv = (rv_l * v_l + p_l) * ny;
      flux_rw = rvw_l * ny;
      flux_rE = ep_l * v_l * ny;
      //
      // A Entropic modification
      //
      p_r = gam1 * (rE_r - 0.5 * (ru_r * ru_r + rv_r * rv_r + rw_r * rw_r) * ir_r);
      lambda_l = v_l * ny - A * sqrt(gam * p_l * ir_l);
      lambda_r = v_r * ny - A * sqrt(gam * p_r * ir_r);
      if ((lambda_l < 0.) && (lambda_r > 0.)) {
        eig_val1 = lambda_l * (lambda_r - eig_val1) / (lambda_r - lambda_l);
      };
      //
      if (eig_val1 < 0.0) {
        //.. In this case A+ is obtained by multiplying the first
        //.. columne of T-1 with the first row of T with eig_val1
        delta_r  = r_r  - r_l;
        delta_ru = ru_r - ru_l;
        delta_rv = rv_r - rv_l;
        delta_rw = rw_r - rw_l;
        delta_rE = rE_r - rE_l;
        //
        scal = scal * iA;
        hny = ny * iA;
        usc = 1.0/c_speed;
        tempo11 = gam1 * usc;
        //.. First colomne of the matrix T-1
        a1 = usc;
        a2 = u_ave * usc;
        a3 = v_ave * usc - hny;
        a4 = w_ave * usc;
        a5 = 0.5 * u2pv2pw2 * usc + 2.5 * c_speed - scal;
        //.. First row of the matrix T * eig_val1
        b1 = 0.5 * (0.5 * tempo11 * u2pv2pw2 + scal);
        b2 = -0.5 * (tempo11 * u_ave);
        b3 = -0.5 * (hny + tempo11 * v_ave);
        b4 = -0.5 * (tempo11 * w_ave);
        b5 = 0.5 * tempo11;
        //
        alpha1 = b1 * delta_r;
        alpha2 = b2 * delta_ru;
        alpha3 = b3 * delta_rv;
        alpha4 = b4 * delta_rw;
        alpha5 = b5 * delta_rE;
        alpha  = eig_val1 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
        //
        flux_r  += a1 * alpha;
        flux_ru += a2 * alpha;
        flux_rv += a3 * alpha;
        flux_rw += a4 * alpha;
        flux_rE += a5 * alpha;
      };
    };

    // Add to interface residuals
    flux[0]  += flux_r;
    flux[1] += flux_ru;
    flux[2] += flux_rv;
    flux[3] += flux_rw;
    flux[4] += flux_rE;
  }

  template <typename PATCH> CUDA_DH void zField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJZ;
    COMPRESSIBLE_RIGHT_PROJZ;

    real alpha1, alpha2, alpha3, alpha4, alpha5, alpha;
    real a1, a2, a3, a4, a5, b1, b2, b3, b4, b5;
    real ep_r, ep_l;
    real ruv_r, ruv_l, ruw_r, ruw_l, rvw_r, rvw_l;
    real lambda_l, lambda_r;
    real delta_r, delta_ru, delta_rv, delta_rw, delta_rE;
    real hnx, hny, hnz;
    real tempo11, usc;

    real flux_r  = 0;
    real flux_ru = 0;
    real flux_rv = 0;
    real flux_rw = 0;
    real flux_rE = 0;

    real nx = 0;
    real ny = 0;
    real nz = A;

    real gam_l = TGas::gamma(var_l);
    real gam_r = TGas::gamma(var_r);
    real gam   = 0.5*(gam_l + gam_r);
    real gam1  = gam - 1.0;

    //
    //.. Compute the ROE Averages
    //
    //.... some useful quantities
    real coef1 = sqrt(r_l);
    real coef2 = sqrt(r_r);
    real somme_coef = coef1 + coef2;
    real isomme_coef = 1.0/somme_coef;
    real h_l = (gam * rE_l*ir_l) - (.5*gam1) * (u_l * u_l + v_l * v_l + w_l * w_l);
    real h_r = (gam * rE_r*ir_r) - (.5*gam1) * (u_r * u_r + v_r * v_r + w_r * w_r);

    //.... averages
    //real rho_ave = coef1 * coef2;
    real u_ave = (coef1 * u_l + coef2 * u_r) * isomme_coef;
    real v_ave = (coef1 * v_l + coef2 * v_r) * isomme_coef;
    real w_ave = (coef1 * w_l + coef2 * w_r) * isomme_coef;
    real h_ave = (coef1 * h_l + coef2 * h_r) * isomme_coef;

    //
    //.. Compute Speed of sound
    real scal = u_ave * nx + v_ave * ny + w_ave * nz;
    real iA = 1.0/A;
    real u2pv2pw2 = u_ave * u_ave + v_ave * v_ave + w_ave * w_ave;
    real c_speed = gam1 * (h_ave - 0.5 * u2pv2pw2);
    if(c_speed < 1e-6) c_speed = 1e-6;    // avoid division by 0 if critical
    c_speed = sqrt(c_speed);
    real c_speed2 = c_speed * A;

    //
    //.. Compute the eigenvalues of the Jacobian matrix
    real eig_val1 = scal - c_speed2;
    real eig_val2 = scal;
    real eig_val3 = scal + c_speed2;
    //
    //.. Compute the ROE flux
    //.... In this part many tests upon the eigenvalues
    //.... are done to simplify calculations
    //.... Here we use the two formes of the ROE flux :
    //.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
    //.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
    //
    if (eig_val2 <= 0.0) {
      p_r = gam1 * (rE_r - 0.5 * (ru_r * ru_r + rv_r * rv_r + rw_r * rw_r) * ir_r);
      ep_r = rE_r + p_r;
      ruv_r = ru_r * v_r;
      ruw_r = ru_r * w_r;
      rvw_r = rv_r * w_r;
      flux_r  = ru_r * nx + rv_r * ny + rw_r * nz;
      flux_ru = (ru_r * u_r + p_r) * nx + ruv_r * ny + ruw_r * nz;
      flux_rv = ruv_r * nx + (rv_r * v_r + p_r) * ny + rvw_r * nz;
      flux_rw = ruw_r * nx + rvw_r * ny + (rw_r * w_r + p_r) * nz;
      flux_rE = ep_r * (u_r * nx + v_r * ny + w_r * nz);
      //
      //.... A Entropic modification
      //
      p_l = gam1 * (rE_l - 0.5 * (ru_l * ru_l + rv_l * rv_l + rw_l * rw_l) * ir_l);
      lambda_l = u_l * nx + v_l * ny + w_l * nz + A * sqrt(gam * p_l * ir_l);
      lambda_r = u_r * nx + v_r * ny + w_r * nz + A * sqrt(gam * p_r * ir_r);
      if ((lambda_l < 0.) && (lambda_r > 0.)) {
        eig_val3 = lambda_r * (eig_val3 - lambda_l) / (lambda_r - lambda_l);
      };
      //
      if (eig_val3 > 0.0) {
        //.. In this case A+ is obtained by multiplying the last
        //.. colomne of T-1 with the last row of T with eig_val3
        delta_r  = r_r  - r_l;
        delta_ru = ru_r - ru_l;
        delta_rv = rv_r - rv_l;
        delta_rw = rw_r - rw_l;
        delta_rE = rE_r - rE_l;
        //
        scal = scal * iA;
        hnx = nx * iA;
        hny = ny * iA;
        hnz = nz * iA;
        usc = 1.0/c_speed;
        tempo11 = gam1 * usc;
        //.. Last columne of the matrix T-1
        a1 = usc;
        a2 = u_ave * usc + hnx;
        a3 = v_ave * usc + hny;
        a4 = w_ave * usc + hnz;
        a5 = 0.5 * u2pv2pw2 * usc + 2.5 * c_speed + scal;
        //.. Last row of the matrix T * eig_val3
        b1 = 0.5 * (0.5 * tempo11 * u2pv2pw2 - scal);
        b2 = 0.5 * (hnx - tempo11 * u_ave);
        b3 = 0.5 * (hny - tempo11 * v_ave);
        b4 = 0.5 * (hnz - tempo11 * w_ave);
        b5 = 0.5 * tempo11;
        //
        alpha1 = b1 * delta_r;
        alpha2 = b2 * delta_ru;
        alpha3 = b3 * delta_rv;
        alpha4 = b4 * delta_rw;
        alpha5 = b5 * delta_rE;
        alpha  = eig_val3 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
        //
        flux_r  -= a1 * alpha;
        flux_ru -= a2 * alpha;
        flux_rv -= a3 * alpha;
        flux_rw -= a4 * alpha;
        flux_rE -= a5 * alpha;
      };
    };
    //
    if(eig_val2 > 0.0) {
      p_l = gam1 * (rE_l - 0.5 * (ru_l * ru_l + rv_l * rv_l + rw_l * rw_l) * ir_l);
      ep_l = rE_l + p_l;
      ruv_l = ru_l * v_l;
      ruw_l = ru_l * w_l;
      rvw_l = rv_l * w_l;
      flux_r  = ru_l * nx + rv_l * ny + rw_l * nz;
      flux_ru = (ru_l * u_l + p_l) * nx + ruv_l * ny + ruw_l * nz;
      flux_rv = ruv_l * nx + (rv_l * v_l + p_l) * ny + rvw_l * nz;
      flux_rw = ruw_l * nx + rvw_l * ny + (rw_l * w_l + p_l) * nz;
      flux_rE = ep_l * (u_l * nx + v_l * ny + w_l * nz);
      //
      // A Entropic modification
      //
      p_r = gam1 * (rE_r - 0.5 * (ru_r * ru_r + rv_r * rv_r + rw_r * rw_r) * ir_r);
      lambda_l = u_l * nx + v_l * ny + w_l * nz - A * sqrt(gam * p_l * ir_l);
      lambda_r   = u_r * nx + v_r * ny + w_r * nz - A * sqrt(gam * p_r * ir_r);
      if ((lambda_l < 0.) && (lambda_r > 0.)) {
        eig_val1 = lambda_l * (lambda_r - eig_val1) / (lambda_r - lambda_l);
      };
      //
      if (eig_val1 < 0.0) {
        //.. In this case A+ is obtained by multiplying the first
        //.. columne of T-1 with the first row of T with eig_val1
        delta_r  = r_r  - r_l;
        delta_ru = ru_r - ru_l;
        delta_rv = rv_r - rv_l;
        delta_rw = rw_r - rw_l;
        delta_rE = rE_r - rE_l;
        //
        scal = scal * iA;
        hnx = nx * iA;
        hny = ny * iA;
        hnz = nz * iA;
        usc = 1.0/c_speed;
        tempo11 = gam1 * usc;
        //.. First colomne of the matrix T-1
        a1 = usc;
        a2 = u_ave * usc - hnx;
        a3 = v_ave * usc - hny;
        a4 = w_ave * usc - hnz;
        a5 = 0.5 * u2pv2pw2 * usc + 2.5 * c_speed - scal;
        //.. First row of the matrix T * eig_val1
        b1 = 0.5 * (0.5 * tempo11 * u2pv2pw2 + scal);
        b2 = -0.5 * (hnx + tempo11 * u_ave);
        b3 = -0.5 * (hny + tempo11 * v_ave);
        b4 = -0.5 * (hnz + tempo11 * w_ave);
        b5 = 0.5 * tempo11;
        //
        alpha1 = b1 * delta_r;
        alpha2 = b2 * delta_ru;
        alpha3 = b3 * delta_rv;
        alpha4 = b4 * delta_rw;
        alpha5 = b5 * delta_rE;
        alpha  = eig_val1 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
        //
        flux_r  += a1 * alpha;
        flux_ru += a2 * alpha;
        flux_rv += a3 * alpha;
        flux_rw += a4 * alpha;
        flux_rE += a5 * alpha;
      };
    };

    // Add to interface residuals
    flux[0]  += flux_r;
    flux[1] += flux_ru;
    flux[2] += flux_rv;
    flux[3] += flux_rw;
    flux[4] += flux_rE;
  }

};

#endif // ROE_H
