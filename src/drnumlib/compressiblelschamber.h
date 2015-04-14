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

template <unsigned int DIM, typename TGas, intreal_t T_press, intreal_t T_temp, intreal_t T_U, intreal_t T_V, intreal_t T_W>
struct CompressibleLsChamber
{

  CUDA_DH static void operate(real* var_f, real*, real* var, real, real, real, real, real, real, real)
  {
    real p_f, T_f, u_f, v_f, w_f;
    TGas::conservativeToPrimitive(var_f, p_f, T_f, u_f, v_f, w_f);

    real p0  = INTREAL(T_press);
    real T0  = INTREAL(T_temp);
    real gam = TGas::gamma(var_f);
    real Ma  = sqrt((pow(p0/p_f, (gam-1)/gam) - 1)*2/(gam-1));
    real p   = p_f;
    real T   = T0/(1 + 0.5*(gam-1)*Ma*Ma);
    real a   = sqrt(gam*TGas::R(var_f)*T);
    real u   = INTREAL(T_U)*a*Ma;
    real v   = INTREAL(T_V)*a*Ma;
    real w   = INTREAL(T_W)*a*Ma;

    TGas::primitiveToConservative(p, T, u, v, w, var);

    if (isnan(var[0]) || isinf(var[0])) {
      printf("Ma=%f, p=%f, T=%f, a=%f, p0=%f, T0=%f, u=%f, v=%f, w=%f, rho=%f\n", Ma, p, T, a, p0, T0, u, v, w, var[0]);
      printf("field cell --> p=%f, T=%f, u=%f, v=%f, w=%f\n\n", p_f, T_f, u_f, v_f, w_f);
    }

  }

};
#endif // COMPRESSIBLELSCHAMBER_H
