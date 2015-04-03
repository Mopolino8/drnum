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

#ifndef COMPRESSIBLELSSLIP_H
#define COMPRESSIBLELSSLIP_H

#include "drnum.h"

template <unsigned int DIM, typename TPatch, typename TGas>
struct CompressibleLsSlip
{
  CUDA_DH static bool usePre() { return false; }
  CUDA_DH static bool usePost() { return false; }

  CUDA_DH static void pre (real* var, real gx, real gy, real gz) {}
  CUDA_DH static void post(real* var, real gx, real gy, real gz)
  {
    real p, T, u, v, w;
    TGas::conservativeToPrimitive(var, p, T, u, v, w);
    real s = u*gx + v*gy + w*gz;
    u -= s*gx;
    v -= s*gy;
    w -= s*gz;
    TGas::primitiveToConservative(p, T, u, v, w, var);
  }

  CUDA_DH static void operate(real* var, real, real wgt, real gx, real gy, real gz)
  {
    real p, T, u, v, w;
    TGas::conservativeToPrimitive(var, p, T, u, v, w);
    real s = u*gx + v*gy + w*gz;
    u -= s*(1 + wgt)*gx;
    v -= s*(1 + wgt)*gy;
    w -= s*(1 + wgt)*gz;
    TGas::primitiveToConservative(p, T, u, v, w, var);
  }

};

#endif // COMPRESSIBLELSSLIP_H
