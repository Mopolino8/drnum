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
#ifndef VANLEER_H
#define VANLEER_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class VanLeer : public CompressibleFlux
{

  TReconstruction* m_Reconstruction;

public: // methods

  VanLeer(TReconstruction* reconstruction) { m_Reconstruction = reconstruction; }

  void xField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);
  void yField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);
  void zField(CartesianPatch *patch,
              size_t i, size_t j, size_t k,
              real x, real y, real z,
              real A, real* flux);


};


template <typename TReconstruction, typename TGas>
void VanLeer<TReconstruction, TGas>::xField(CartesianPatch *patch,
                                            size_t i, size_t j, size_t k,
                                            real x, real y, real z,
                                            real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJX;
  COMPRESSIBLE_RIGHT_PROJX;

  real M_l = u_l/a_l;
  real M_r = u_r/a_r;
  countFlops(2);

  real F0_l = 0, F1_l = 0, F2_l = 0, F3_l = 0, F4_l = 0;
  if (M_l >= 1) {
    F0_l = r_l*u_l;
    F1_l = u_l*ru_l + p_l;
    F2_l = u_l*rv_l;
    F3_l = u_l*rw_l;
    F4_l = u_l*(rE_l + p_l);
    countFlops(7);
  } else if (M_l > -1) {
    F0_l = 0.25*a_l*r_l*(M_l+1)*(M_l+1);
    F1_l = F0_l*(u_l + p_l/(a_l*r_l)*(-M_l+2));
    F2_l = F0_l*v_l;
    F3_l = F0_l*w_l;
    F4_l = F0_l/r_l*(rE_l + p_l);
    countFlops(18);
  };

  real F0_r = 0, F1_r = 0, F2_r = 0, F3_r = 0, F4_r = 0;
  if (M_r <= -1) {
    F0_r = r_r*u_r;
    F1_r = u_r*ru_r + p_r;
    F2_r = u_r*rv_r;
    F3_r = u_r*rw_r;
    F4_r = u_r*(rE_r + p_r);
    countFlops(7);
  } else if (M_r < 1) {
    F0_r = -0.25*a_r*r_r*(M_r-1)*(M_r-1);
    F1_r = F0_r*(u_r + p_r/(a_r*r_r)*(-M_r-2));
    F2_r = F0_r*v_r;
    F3_r = F0_r*w_r;
    F4_r = F0_r/r_r*(rE_r + p_r);
    countFlops(18);
  };

  flux[0] += A*(F0_r + F0_l);
  flux[1] += A*(F1_r + F1_l);
  flux[2] += A*(F2_r + F2_l);
  flux[3] += A*(F3_r + F3_l);
  flux[4] += A*(F4_r + F4_l);
  countFlops(10);
}

template <typename TReconstruction, typename TGas>
void VanLeer<TReconstruction, TGas>::yField(CartesianPatch *patch,
                                            size_t i, size_t j, size_t k,
                                            real x, real y, real z,
                                            real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJY;
  COMPRESSIBLE_RIGHT_PROJY;

  real M_l = v_l/a_l;
  real M_r = v_r/a_r;
  countFlops(2);

  real F0_l = 0, F1_l = 0, F2_l = 0, F3_l = 0, F4_l = 0;
  if (M_l >= 1) {
    F0_l = r_l*v_l;
    F1_l = v_l*ru_l;
    F2_l = v_l*rv_l + p_l;
    F3_l = v_l*rw_l;
    F4_l = v_l*(rE_l + p_l);
    countFlops(7);
  } else if (M_l > -1) {
    F0_l = 0.25*a_l*r_l*(M_l+1)*(M_l+1);
    F1_l = F0_l*u_l;
    F2_l = F0_l*(v_l + p_l/(a_l*r_l)*(-M_l+2));
    F3_l = F0_l*w_l;
    F4_l = F0_l/r_l*(rE_l + p_l);
    countFlops(18);
  };

  real F0_r = 0, F1_r = 0, F2_r = 0, F3_r = 0, F4_r = 0;
  if (M_r <= -1) {
    F0_r = r_r*v_r;
    F1_r = v_r*ru_r;
    F2_r = v_r*rv_r + p_r;
    F3_r = v_r*rw_r;
    F4_r = v_r*(rE_r + p_r);
    countFlops(7);
  } else if (M_r < 1) {
    F0_r = -0.25*a_r*r_r*(M_r-1)*(M_r-1);
    F1_r = F0_r*u_r;
    F2_r = F0_r*(v_r + p_r/(a_r*r_r)*(-M_r-2));
    F3_r = F0_r*w_r;
    F4_r = F0_r/r_r*(rE_r + p_r);
    countFlops(18);
  };

  flux[0] += A*(F0_r + F0_l);
  flux[1] += A*(F1_r + F1_l);
  flux[2] += A*(F2_r + F2_l);
  flux[3] += A*(F3_r + F3_l);
  flux[4] += A*(F4_r + F4_l);
  countFlops(10);
}

template <typename TReconstruction, typename TGas>
void VanLeer<TReconstruction, TGas>::zField(CartesianPatch *patch,
                                            size_t i, size_t j, size_t k,
                                            real x, real y, real z,
                                            real A, real* flux)
{
  COMPRESSIBLE_LEFT_PROJZ;
  COMPRESSIBLE_RIGHT_PROJZ;

  real M_l = w_l/a_l;
  real M_r = w_r/a_r;
  countFlops(2);

  real F0_l = 0, F1_l = 0, F2_l = 0, F3_l = 0, F4_l = 0;
  if (M_l >= 1) {
    F0_l = r_l*w_l;
    F1_l = w_l*ru_l;
    F2_l = w_l*rv_l;
    F3_l = w_l*rw_l + p_l;
    F4_l = w_l*(rE_l + p_l);
    countFlops(7);
  } else if (M_l > -1) {
    F0_l = 0.25*a_l*r_l*(M_l+1)*(M_l+1);
    F1_l = F0_l*u_l;
    F2_l = F0_l*v_l;
    F3_l = F0_l*(w_l + p_l/(a_l*r_l)*(-M_l+2));
    F4_l = F0_l/r_l*(rE_l + p_l);
    countFlops(18);
  };

  real F0_r = 0, F1_r = 0, F2_r = 0, F3_r = 0, F4_r = 0;
  if (M_r <= -1) {
    F0_r = r_r*w_r;
    F1_r = w_r*ru_r;
    F2_r = w_r*rv_r;
    F3_r = w_r*rw_r + p_r;
    F4_r = w_r*(rE_r + p_r);
    countFlops(7);
  } else if (M_r < 1) {
    F0_r = -0.25*a_r*r_r*(M_r-1)*(M_r-1);
    F1_r = F0_r*u_r;
    F2_r = F0_r*v_r;
    F3_r = F0_r*(w_r + p_r/(a_r*r_r)*(-M_r-2));
    F4_r = F0_r/r_r*(rE_r + p_r);
    countFlops(18);
  };

  flux[0] += A*(F0_r + F0_l);
  flux[1] += A*(F1_r + F1_l);
  flux[2] += A*(F2_r + F2_l);
  flux[3] += A*(F3_r + F3_l);
  flux[4] += A*(F4_r + F4_l);
  countFlops(10);
}

#endif // VANLEER_H
