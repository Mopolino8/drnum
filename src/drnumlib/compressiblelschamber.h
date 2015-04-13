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

#ifndef COMPRESSIBLELSCHAMBER_H
#define COMPRESSIBLELSCHAMBER_H

#include "drnum.h"

template <unsigned int DIM, typename TPatch, typename TGas, int T_temp, int T_press>
struct CompressibleLsChamber
{
  CUDA_DH static real extrapolate(real h0, real h1, real h2, real f1, real f2)
  {
    return f1 + (h0-h1)*(f1-f2)/(h1-h2);
  }

  CUDA_DH static void operate(real* var_f, real*, real* var, real h0, real h1, real h2, real wgt, real gx, real gy, real gz)
  {
    real p_f, T_f, u_f, v_f, w_f;
    real p_bc, T_bc, u_bc, v_bc, w_bc;
    TGas::conservativeToPrimitive(var_f, p_f, T_f, u_f, v_f, w_f);

    //real a_f = CHECKED_REAL(sqrt(TGas::gamma(var_f)*TGas::R(var_f)*T_f));
    //real


    TGas::primitiveToConservative(p, T, u, v, w, var);
  }

};
#endif // COMPRESSIBLELSCHAMBER_H
