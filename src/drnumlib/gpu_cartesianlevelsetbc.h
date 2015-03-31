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

#ifndef GPU_CARTESIANLEVELSETBC_H
#define GPU_CARTESIANLEVELSETBC_H

template <unsigned int DIM, typename BC>
class GPU_CartesianLevelSetBC;

#include "gpu_levelsetbc.h"

#include "drnum.h"
#include "genericoperation.h"
#include "gpu_cartesianpatch.h"

template <unsigned int DIM, unsigned int IVAR, typename BC>
class GPU_CartesianLevelSetBC : public GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch, BC>
{

protected: // attributes


public: // methods

  GPU_LevelSetBC(PatchGrid* patch_grid);

  CUDA_DO static real  g     (GPU_CartesianPatch &patch, size_t i, size_t j, size_t k);
  CUDA_DO static real& marker(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k);
  CUDA_DO static bool  marked(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k);
  CUDA_DO static void  grad  (GPU_CartesianPatch &patch, size_t i, size_t j, size_t k, real &gx, real &gy, real &gz);

};

//CUDA_DH real& f(size_t i_field, size_t i_var, size_t i, size_t j, size_t k)


template <unsigned int DIM, typename BC>
real GPU_CartesianLevelSetBC<DIM,BC>::g(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k)
{
  return patch.f(0, IVAR, i, j, k);
}

template <unsigned int DIM, typename BC>
real& GPU_CartesianLevelSetBC<DIM,BC>::marker(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k)
{
  return patch.f(2, 0, i, j, k);
}

template <unsigned int DIM, typename BC>
bool GPU_CartesianLevelSetBC<DIM,BC>::marked(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k)
{
  return marker(patch, i, j, k) > 0.5;
}

template <unsigned int DIM, typename BC>
void GPU_CartesianLevelSetBC<DIM,BC>::grad(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k, real &gx, real &gy, real &gz)
{
  gx = 0.5*patch.idx()*(g(patch, i+1, j, k) - g(patch, i-1, j, k));
  gy = 0.5*patch.idy()*(g(patch, i, j+1, k) - g(patch, i, j+1, k));
  gz = 0.5*patch.idz()*(g(patch, i, j, k+1) - g(patch, i, j, k+1));
}

template <unsigned int DIM, typename BC>
GPU_CartesianLevelSetBC<DIM,BC>::GPU_CartesianLevelSetBC(PatchGrid* patch_grid)
  : GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch, BC>(patch_grid)
{
}


template <unsigned int DIM, typename BC>
__global__ GPU_CartesianLevelSetBC_Operate(GPU_CartesianPatch patch)
{
  size_t i = 1 + blockIdx.x;
  size_t j = 1 + blockIdx.y;
  size_t k = 1 + threadIdx.x;

  if (i >= (patch.sizeI() - 1) || j >= (patch.sizeJ() - 1) || k >= (patch.sizeK() - 1)) {
    return;
  }


  if (g(i, j, k) < 0) {
    real gx_c, gy_c, gz_c;
    GPU_CartesianLevelSetBC<DIM,BC>::grad(patch, i, j, k, gx_c, gy_c, gz_c);
    if (g(i+1, j, k) >= 0) {
      real gx = patch.idx()*(GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i+1, j, k) - GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k));
      m_Bc.operate(i, j, k, i+1, j, k, gx, gy_c, gz_c);
    }
    if (g(i-1, j, k) >= 0) {
      real gx = patch.idx()*(GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k) - GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i-1, j, k));
      m_Bc.operate(i, j, k, i-1, j, k, gx, gy_c, gz_c);
    }
    if (g(i, j+1, k) >= 0) {
      real gy = patch.idy()*(GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j+1, k) - GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k));
      m_Bc.operate(i, j, k, i, j+1, k, gx_c, gy, gz_c);
    }
    if (g(i, j-1, k) >= 0) {
      real gy = patch.idy()*(GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k) - GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j-1, k));
      m_Bc.operate(i, j, k, i, j-1, k, gx_c, gy, gz_c);
    }
    if (g(i, j, k+1) >= 0) {
      real gz = patch.idz()*(GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k+1) - GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k));
      m_Bc.operate(i, j, k, i, j, k+1, gx_c, gy, gz_c);
    }
    if (g(i, j, k-1) >= 0) {
      real gz = patch.idz()*(GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k) - GPU_CartesianLevelSetBC<DIM,BC>::g(patch, i, j, k-1));
      m_Bc.operate(i, j, k, i, j, k-1, gx_c, gy, gz_c);
    }
  }
}

#endif // GPU_CARTESIANLEVELSETBC_H
