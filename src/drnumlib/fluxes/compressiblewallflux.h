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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef COMPRESSIBLEWALLFLUX_H
#define COMPRESSIBLEWALLFLUX_H

#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class CompressibleWallFlux
{

  TReconstruction m_Reconstruction;

public:

  template <typename PATCH> CUDA_DH void xWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i-1, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*p;
    flux[2] += A*mu*v*patch->idx();
    flux[3] += A*mu*w*patch->idx();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void yWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j-1, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*mu*u*patch->idy();
    flux[2] += A*p;
    flux[3] += A*mu*w*patch->idy();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void zWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k-1, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*mu*u*patch->idz();
    flux[2] += A*mu*v*patch->idz();
    flux[3] += A*p;
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void xWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*p;
    flux[2] -= A*mu*v*patch->idx();
    flux[3] -= A*mu*w*patch->idx();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void yWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] -= A*mu*u*patch->idy();
    flux[2] += A*p;
    flux[3] -= A*mu*w*patch->idy();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void zWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] -= A*mu*u*patch->idz();
    flux[2] -= A*mu*v*patch->idz();
    flux[3] += A*p;
    countFlops(2);
  }

};

#endif // COMPRESSIBLEWALLFLUX_H
