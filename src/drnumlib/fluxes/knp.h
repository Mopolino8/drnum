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
#ifndef KNP_H
#define KNP_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class KNP : public CompressibleFlux
{

  TReconstruction m_Reconstruction;

public: // methods

  template <typename PATCH> CUDA_DH void xField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJX;
    COMPRESSIBLE_RIGHT_PROJX;

    real psi_l = max(max(a_l + u_l, a_r + u_r), real(0.0));
    real psi_r = max(max(a_l - u_l, a_r - u_r), real(0.0));
    real alpha = psi_l/(psi_l + psi_r);
    real omega = alpha*(1 - alpha)*(psi_l + psi_r);
    countFlops(10);

    real F0_l = r_l;
    real F1_l = ru_l;
    real F2_l = rv_l;
    real F3_l = rw_l;
    real F4_l = rE_l + p_l;
    countFlops(1);

    real F0_r = r_r;
    real F1_r = ru_r;
    real F2_r = rv_r;
    real F3_r = rw_r;
    real F4_r = rE_r + p_r;
    countFlops(1);

    flux[0] += A*(alpha*u_l*F0_l + (1-alpha)*u_r*F0_r - omega*(F0_r - F0_l));
    flux[1] += A*(alpha*u_l*F1_l + (1-alpha)*u_r*F1_r - omega*(F1_r - F1_l));
    flux[2] += A*(alpha*u_l*F2_l + (1-alpha)*u_r*F2_r - omega*(F2_r - F2_l));
    flux[3] += A*(alpha*u_l*F3_l + (1-alpha)*u_r*F3_r - omega*(F3_r - F3_l));
    flux[4] += A*(alpha*u_l*F4_l + (1-alpha)*u_r*F4_r - omega*(F4_r - F4_l));
    countFlops(55);

    flux[1] += A*(alpha*p_l + (1-alpha)*p_r);
    countFlops(6);
  }

  template <typename PATCH> CUDA_DH void yField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJY;
    COMPRESSIBLE_RIGHT_PROJY;

    real psi_l = max(max(a_l + v_l, a_r + v_r), real(0.0));
    real psi_r = max(max(a_l - v_l, a_r - v_r), real(0.0));
    real alpha = psi_l/(psi_l + psi_r);
    real omega = alpha*(1 - alpha)*(psi_l + psi_r);
    countFlops(10);

    real F0_l = r_l;
    real F1_l = ru_l;
    real F2_l = rv_l;
    real F3_l = rw_l;
    real F4_l = rE_l + p_l;
    countFlops(1);

    real F0_r = r_r;
    real F1_r = ru_r;
    real F2_r = rv_r;
    real F3_r = rw_r;
    real F4_r = rE_r + p_r;
    countFlops(1);

    flux[0] += A*(alpha*v_l*F0_l + (1-alpha)*v_r*F0_r - omega*(F0_r - F0_l));
    flux[1] += A*(alpha*v_l*F1_l + (1-alpha)*v_r*F1_r - omega*(F1_r - F1_l));
    flux[2] += A*(alpha*v_l*F2_l + (1-alpha)*v_r*F2_r - omega*(F2_r - F2_l));
    flux[3] += A*(alpha*v_l*F3_l + (1-alpha)*v_r*F3_r - omega*(F3_r - F3_l));
    flux[4] += A*(alpha*v_l*F4_l + (1-alpha)*v_r*F4_r - omega*(F4_r - F4_l));
    countFlops(55);

    flux[2] += A*(alpha*p_l + (1-alpha)*p_r);
    countFlops(6);
  }

  template <typename PATCH> CUDA_DH void zField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJZ;
    COMPRESSIBLE_RIGHT_PROJZ;

    real psi_l = max(max(a_l + w_l, a_r + w_r), real(0.0));
    real psi_r = max(max(a_l - w_l, a_r - w_r), real(0.0));
    real alpha = psi_l/(psi_l + psi_r);
    real omega = alpha*(1 - alpha)*(psi_l + psi_r);
    countFlops(10);

    real F0_l = r_l;
    real F1_l = ru_l;
    real F2_l = rv_l;
    real F3_l = rw_l;
    real F4_l = rE_l + p_l;
    countFlops(1);

    real F0_r = r_r;
    real F1_r = ru_r;
    real F2_r = rv_r;
    real F3_r = rw_r;
    real F4_r = rE_r + p_r;
    countFlops(1);

    flux[0] += A*(alpha*w_l*F0_l + (1-alpha)*w_r*F0_r - omega*(F0_r - F0_l));
    flux[1] += A*(alpha*w_l*F1_l + (1-alpha)*w_r*F1_r - omega*(F1_r - F1_l));
    flux[2] += A*(alpha*w_l*F2_l + (1-alpha)*w_r*F2_r - omega*(F2_r - F2_l));
    flux[3] += A*(alpha*w_l*F3_l + (1-alpha)*w_r*F3_r - omega*(F3_r - F3_l));
    flux[4] += A*(alpha*w_l*F4_l + (1-alpha)*w_r*F4_r - omega*(F4_r - F4_l));
    countFlops(55);

    flux[3] += A*(alpha*p_l + (1-alpha)*p_r);
    countFlops(6);
  }


};

#endif // KNP_H
