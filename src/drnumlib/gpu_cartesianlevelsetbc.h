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

template <unsigned int DIM, typename LS, typename BC>
class GPU_CartesianLevelSetBC;

#include "gpu_levelsetbc.h"

#include "drnum.h"
#include "genericoperation.h"
#include "gpu_cartesianpatch.h"

template <unsigned int DIM, typename LS, typename BC>
class GPU_CartesianLevelSetBC : public GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>
{

protected: // attributes


public: // methods

  GPU_CartesianLevelSetBC(PatchGrid* patch_grid, int cuda_device = 0, size_t thread_limit = 0);

  CUDA_DO static void grad(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k, real &gx, real &gy, real &gz);
  CUDA_DO static void getOutsideState(GPU_CartesianPatch &patch,
                                      size_t i, size_t j, size_t k,
                                      int di, int dj, int dk,
                                      real gx, real gy, real gz,
                                      real &h, real &w, real* var);

  CUDA_HO virtual void operator()();

};


template <unsigned int DIM, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,LS,BC>::grad(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k, real &gx, real &gy, real &gz)
{
  gx = 0.5*patch.idx()*(LS::G(patch, i+1, j, k) - LS::G(patch, i-1, j, k));
  gy = 0.5*patch.idy()*(LS::G(patch, i, j+1, k) - LS::G(patch, i, j-1, k));
  gz = 0.5*patch.idz()*(LS::G(patch, i, j, k+1) - LS::G(patch, i, j, k-1));
}

template <unsigned int DIM, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,LS,BC>::getOutsideState(GPU_CartesianPatch &patch,
                                                         size_t i, size_t j, size_t k,
                                                         int di, int dj, int dk,
                                                         real gx, real gy, real gz,
                                                         real &h, real &w, real *var)
{
  
  // careful with parallel level sets (e.g. flat plate)
  // there might be a 0/0 occurring
  
  dim_t<DIM> dim;
  real h0 = LS::G(patch, i, j, k);
  real h1 = LS::G(patch, i + di, j + dj, k + dk);
  real h2 = LS::G(patch, i + 2*di, j + 2*dj, k + 2*dk);
  if (h2 < 0) {
    patch.getVar(dim, 0, i + di, j + dj, k + dk, var);
    h = h1;
    w = 0;
  } else {
    patch.getVar(dim, 0, i + di, j + dj, k + dk, var);
    real var2[DIM];
    patch.getVar(dim, 0, i + 2*di, j + 2*dj, k + 2*dk, var2);
    w = h1/(h1 - h0);
    h = w*h1 + (1 - w)*h2;
    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = w*var[i_var] + (1 - w)*var2[i_var];
    }
    w = -h0/h;
  }
}

template <unsigned int DIM, typename LS, typename BC>
GPU_CartesianLevelSetBC<DIM,LS,BC>::GPU_CartesianLevelSetBC(PatchGrid* patch_grid, int cuda_device, size_t thread_limit)
  : GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>(patch_grid, cuda_device, thread_limit)
{
}

#define GPU_CARTESIANLEVELSETBC_INDICES               \
  size_t i = blockIdx.x + 2;                          \
  size_t j = blockDim.y*blockIdx.y + threadIdx.y + 2; \
  size_t k = threadIdx.x + 2;


template <unsigned int DIM, typename LS, typename BC>
__global__ void GPU_CartesianLevelSetBC_kernelOperate(GPU_CartesianPatch patch)
{
  GPU_CARTESIANLEVELSETBC_INDICES;

  dim_t<DIM> dim;

  if (i >= (patch.sizeI() - 2) || j >= (patch.sizeJ() - 2) || k >= (patch.sizeK() - 2)) {
    return;
  }

  if (LS::G(patch, i, j, k) < 0) {

    real var_outside[DIM], var[DIM];
    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = 0.0;
    }

    real gx, gy, gz;
    int count = 0;
    real total_weight = 0;
    real h, w;
    GPU_CartesianLevelSetBC<DIM,LS,BC>::grad(patch, i, j, k, gx, gy, gz);
    for (int di = -1; di <= 1; ++di) {
      for (int dj = -1; dj <= 1; ++dj) {
        for (int dk = -1; dk <= 1; ++dk) {
          if (di != 0 || dj != 0 || dk != 0) {
            if (LS::G(patch, i + di, j + dj, k + dk) >= 0) {
              real dx = di*patch.dx();
              real dy = dj*patch.dy();
              real dz = dk*patch.dz();
              real weight = fabs(dx*gx + dy*gy + dz*gz)/sqrt(dx*dx + dy*dy + dz*dz);
              GPU_CartesianLevelSetBC<DIM,LS,BC>::getOutsideState(patch, i, j, k, di, dj, dk, gx, gy, gz, h, w, var_outside);
              BC::operate(var_outside, h, w, gx, gy, gz);
              total_weight += weight;
              ++count;
              for (size_t i_var = 0; i_var < DIM; ++i_var) {
                var[i_var] += weight*var_outside[i_var];
              }
            }
          }
        }
      }
    }
    if (count > 0) {
      for (size_t i_var = 0; i_var < DIM; ++i_var) {
        var[i_var] /= total_weight;
      }
      real G_self = LS::G(patch, i, j, k);
      patch.setVar(dim, 0, i, j, k, var);
      LS::updateG(patch, i, j, k, G_self);
    }

  }
}

template <unsigned int DIM, typename LS, typename BC>
__global__ void GPU_CartesianLevelSetBC_kernelExtrapolate(GPU_CartesianPatch patch)
{
  GPU_CARTESIANLEVELSETBC_INDICES;

  dim_t<DIM> dim;

  if (i >= (patch.sizeI() - 2) || j >= (patch.sizeJ() - 2) || k >= (patch.sizeK() - 2)) {
    return;
  }

  real G_self = LS::G(patch, i, j, k);

  if (G_self < 0) {

    real var_neigh[DIM], var[DIM];
    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = 0.0;
    }

    real gx, gy, gz;
    int count1 = 0;
    int count2 = 0;
    real total_weight = 0;
    GPU_CartesianLevelSetBC<DIM,LS,BC>::grad(patch, i, j, k, gx, gy, gz);
    for (int di = -1; di <= 1; ++di) {
      for (int dj = -1; dj <= 1; ++dj) {
        for (int dk = -1; dk <= 1; ++dk) {
          if (di != 0 || dj != 0 || dk != 0) {
            if (LS::G(patch, i + di, j + dj, k + dk) < 0) {
              if (LS::G(patch, i + di, j + dj, k + dk) > G_self) {
                real dx = di*patch.dx();
                real dy = dj*patch.dy();
                real dz = dk*patch.dz();
                real weight = fabs(dx*gx + dy*gy + dz*gz)/sqrt(dx*dx + dy*dy + dz*dz);
                total_weight += weight;
                patch.getVar(dim, 0, i + di, j + dj, k + dk, var_neigh);
                ++count1;
                for (size_t i_var = 0; i_var < DIM; ++i_var) {
                  var[i_var] += weight*var_neigh[i_var];
                }
              }
            } else {
              ++count2;
            }
          }
        }
      }
    }
    if (count1 > 0 && count2 == 0) {
      for (size_t i_var = 0; i_var < DIM; ++i_var) {
        var[i_var] /= total_weight;
      }
      real G_self = LS::G(patch, i, j, k);
      patch.setVar(dim, 2, i, j, k, var);
      LS::updateG(patch, i, j, k, G_self, 2);
    }

  }
}

template <unsigned int DIM, typename LS, typename BC>
__global__ void GPU_CartesianLevelSetBC_kernelPre(GPU_CartesianPatch patch)
{
  GPU_CARTESIANLEVELSETBC_INDICES;

  dim_t<DIM> dim;

  if (i >= (patch.sizeI() - 1) || j >= (patch.sizeJ() - 1) || k >= (patch.sizeK() - 1)) {
    return;
  }

  if (LS::G(patch, i, j, k) >= 0) {

    real var[DIM];
    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = 0.0;
    }

    real gx, gy, gz;
    bool crossover = false;
    GPU_CartesianLevelSetBC<DIM,LS,BC>::grad(patch, i, j, k, gx, gy, gz);
    if      (LS::G(patch, i+1, j, k) < 0) crossover = true;
    else if (LS::G(patch, i-1, j, k) < 0) crossover = true;
    else if (LS::G(patch, i, j+1, k) < 0) crossover = true;
    else if (LS::G(patch, i, j-1, k) < 0) crossover = true;
    else if (LS::G(patch, i, j, k+1) < 0) crossover = true;
    else if (LS::G(patch, i, j, k-1) < 0) crossover = true;

    if (crossover) {
      BC::pre(var, gx, gy, gz);
      patch.setVar(dim, 0, i, j, k, var);
    }

  }
}

template <unsigned int DIM, typename LS, typename BC>
__global__ void GPU_CartesianLevelSetBC_kernelPost(GPU_CartesianPatch patch)
{
  GPU_CARTESIANLEVELSETBC_INDICES;

  dim_t<DIM> dim;

  if (i >= (patch.sizeI() - 1) || j >= (patch.sizeJ() - 1) || k >= (patch.sizeK() - 1)) {
    return;
  }

  if (LS::G(patch, i, j, k) >= 0) {

    real var[DIM];
    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = 0.0;
    }

    real gx, gy, gz;
    bool crossover = false;
    if      (LS::G(patch, i+1, j, k) < 0) crossover = true;
    else if (LS::G(patch, i-1, j, k) < 0) crossover = true;
    else if (LS::G(patch, i, j+1, k) < 0) crossover = true;
    else if (LS::G(patch, i, j-1, k) < 0) crossover = true;
    else if (LS::G(patch, i, j, k+1) < 0) crossover = true;
    else if (LS::G(patch, i, j, k-1) < 0) crossover = true;

    if (crossover) {
      GPU_CartesianLevelSetBC<DIM,LS,BC>::grad(patch, i, j, k, gx, gy, gz);
      BC::post(var, gx, gy, gz);
      patch.setVar(dim, 0, i, j, k, var);
    }

  }
}

template <unsigned int DIM, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,LS,BC>::operator()()
{
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  CUDA_CHECK_ERROR;

  //GPU_LevelSetBC<DIM,CartesianPatch,GPU_CartesianPatch,BC>::copyField(0, 2);
  //this->copyField(0, 2);
  cudaDeviceSynchronize();

  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {

    CUDA_CHECK_ERROR;

    size_t max_num_threads = this->m_MaxNumThreads;
    size_t k_lines = max(size_t(1), size_t(max_num_threads/this->m_Patches[i_patch]->sizeK()));

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), this->m_Patches[i_patch]->sizeJ()/k_lines+1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), k_lines, 1);

      if (BC::usePre()) {
        GPU_CartesianLevelSetBC_kernelPre<DIM,LS,BC> <<<blocks, threads>>>(this->m_GpuPatches[i_patch]);
        CUDA_CHECK_ERROR;
        cudaDeviceSynchronize();
      }

      GPU_CartesianLevelSetBC_kernelOperate<DIM,LS,BC> <<<blocks, threads>>>(this->m_GpuPatches[i_patch]);
      CUDA_CHECK_ERROR;
      cudaDeviceSynchronize();

      if (BC::usePost()) {
        GPU_CartesianLevelSetBC_kernelPost<DIM,LS,BC> <<<blocks, threads>>>(this->m_GpuPatches[i_patch]);
        CUDA_CHECK_ERROR;
        cudaDeviceSynchronize();
      }

      cudaMemcpy(this->m_GpuPatches[i_patch].getField(2), this->m_GpuPatches[i_patch].getField(0), this->m_GpuPatches[i_patch].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
      GPU_CartesianLevelSetBC_kernelExtrapolate<DIM,LS,BC> <<<blocks, threads>>>(this->m_GpuPatches[i_patch]);
      CUDA_CHECK_ERROR;
      cudaDeviceSynchronize();
      cudaMemcpy(this->m_GpuPatches[i_patch].getField(0), this->m_GpuPatches[i_patch].getField(2), this->m_GpuPatches[i_patch].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);

    }

  }
}

#endif // GPU_CARTESIANLEVELSETBC_H
