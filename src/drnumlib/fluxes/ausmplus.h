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
#ifndef AUSMPLUS_H
#define AUSMPLUS_H

#include "fluxes/compressibleflux.h"
#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class AusmPlus : public CompressibleFlux
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

    real a    = FR12*(a_l + a_r);
    real M    = M4(u_l/a, 1) + M4(u_r/a, -1);
    real Mp   = FR12*(M + fabs(M));
    real Mm   = FR12*(M - fabs(M));
    real p    = P5(u_l/a, 1)*p_l + P5(u_r/a, -1)*p_r;
    countFlops(14);

    flux[0] += a*A*(r_l*Mp + r_r*Mm);
    flux[1] += FR12*flux[0]*(u_l + u_r) + A*p - FR12*fabs(flux[0])*(u_r - u_l);
    flux[2] += FR12*flux[0]*(v_l + v_r)       - FR12*fabs(flux[0])*(v_r - v_l);
    flux[3] += FR12*flux[0]*(w_l + w_r)       - FR12*fabs(flux[0])*(w_r - w_l);
    flux[4] += FR12*flux[0]*(H_l + H_r)       - FR12*fabs(flux[0])*(H_r - H_l);
    countFlops(36);
  }

  template <typename PATCH> CUDA_DH void yField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJY;
    COMPRESSIBLE_RIGHT_PROJY;

    real a    = FR12*(a_l + a_r);
    real M    = M4(v_l/a, 1) + M4(v_r/a, -1);
    real Mp   = FR12*(M + fabs(M));
    real Mm   = FR12*(M - fabs(M));
    real p    = P5(v_l/a, 1)*p_l + P5(v_r/a, -1)*p_r;
    countFlops(14);

    flux[0] += a*A*(r_l*Mp + r_r*Mm);
    flux[1] += FR12*flux[0]*(u_l + u_r)       - FR12*fabs(flux[0])*(u_r - u_l);
    flux[2] += FR12*flux[0]*(v_l + v_r) + A*p - FR12*fabs(flux[0])*(v_r - v_l);
    flux[3] += FR12*flux[0]*(w_l + w_r)       - FR12*fabs(flux[0])*(w_r - w_l);
    flux[4] += FR12*flux[0]*(H_l + H_r)       - FR12*fabs(flux[0])*(H_r - H_l);
    countFlops(36);
  }

  template <typename PATCH> CUDA_DH void zField(PATCH *patch,
                                                size_t i, size_t j, size_t k,
                                                real x, real y, real z,
                                                real A, real* flux)
  {
    COMPRESSIBLE_LEFT_PROJZ;
    COMPRESSIBLE_RIGHT_PROJZ;

    real a    = FR12*(a_l + a_r);
    real M    = M4(w_l/a, 1) + M4(w_r/a, -1);
    real Mp   = FR12*(M + fabs(M));
    real Mm   = FR12*(M - fabs(M));
    real p    = P5(w_l/a, 1)*p_l + P5(w_r/a, -1)*p_r;
    countFlops(14);

    flux[0] += a*A*(r_l*Mp + r_r*Mm);
    flux[1] += FR12*flux[0]*(u_l + u_r)       - FR12*fabs(flux[0])*(u_r - u_l);
    flux[2] += FR12*flux[0]*(v_l + v_r)       - FR12*fabs(flux[0])*(v_r - v_l);
    flux[3] += FR12*flux[0]*(w_l + w_r) + A*p - FR12*fabs(flux[0])*(w_r - w_l);
    flux[4] += FR12*flux[0]*(H_l + H_r)       - FR12*fabs(flux[0])*(H_r - H_l);
    countFlops(36);
  }


};


#endif // AUSMPLUS_H
