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
  CUDA_DH static bool usePost() { return true; }

  CUDA_DH static void pre(TPatch* patch, real* var, size_t i, size_t j, size_t k, real gx, real gy, real gz) {}

  CUDA_DH static void operateXYZ(TPatch* patch, real* var, size_t i, size_t j, size_t k, real gx, real gy, real gz)
  {
    dim_t<DIM> dim;
    patch->getVar(dim, 0, i, j, k, var);
    real p, T, u, v, w;
    TGas::conservativeToPrimitive(var, p, T, u, v, w);
    real s = u*gx + v*gy + w*gz;
    u -= s*gx;
    v -= s*gy;
    w -= s*gz;
    TGas::primitiveToConservative(p, T, u, v, w, var);
  }

  CUDA_DH static void post(TPatch* patch, real* var, size_t i, size_t j, size_t k, real gx, real gy, real gz)
  {
    operateXYZ(patch, var, i, j, k, gx, gy, gz);
  }

  CUDA_DH static void operateX(TPatch* patch, real* var, real, size_t i, size_t j, size_t k, size_t dir, real gx, real gy, real gz)
  {
    operateXYZ(patch, var, i+dir, j, k, gx, gy, gz);
  }

  CUDA_DH static void operateY(TPatch* patch, real* var, real, size_t i, size_t j, size_t k, size_t dir, real gx, real gy, real gz)
  {
    operateXYZ(patch, var, i, j+dir, k, gx, gy, gz);
  }

  CUDA_DH static void operateZ(TPatch* patch, real* var, real, size_t i, size_t j, size_t k, size_t dir, real gx, real gy, real gz)
  {
    operateXYZ(patch, var, i, j, k+dir, gx, gy, gz);
  }

};

#endif // COMPRESSIBLELSSLIP_H
