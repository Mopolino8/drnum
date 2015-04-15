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

template <typename TGas>
class CompressibleLsChamber
{

protected: // attributes

  real m_P0;
  real m_T0;


public: // methods

  CUDA_DH CompressibleLsChamber() {}
  CUDA_DH CompressibleLsChamber(real p0, real T0)
  {
    m_P0 = p0;
    m_T0 = T0;
  }

  CUDA_DH void operate(real* var_f, real*, real* var, real, real, real, real, real, real, real)
  {
    real p_f, T_f, u_f, v_f, w_f;
    TGas::conservativeToPrimitive(var_f, p_f, T_f, u_f, v_f, w_f);

    real gam = TGas::gamma(var_f);
    real Ma  = sqrt((pow(m_P0/p_f, (gam-1)/gam) - 1)*2/(gam-1));
    real p   = p_f;
    real T   = m_T0/(1 + 0.5*(gam-1)*Ma*Ma);
    real a   = sqrt(gam*TGas::R(var_f)*T);

    TGas::primitiveToConservative(p, T, Ma*a, 0, 0, var);

  }

};
#endif // COMPRESSIBLELSCHAMBER_H
