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

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
class GPU_CartesianLevelSetBC;

#include "gpu_levelsetbc.h"

#include "drnum.h"
#include "genericoperation.h"
#include "gpu_cartesianpatch.h"

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
class GPU_CartesianLevelSetBC : public GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>
{

public: // data types

  struct cell_t
  {
    size_t i, j, k;
    bool   extrapolate;
  };


protected: // attributes

  QVector<QList<cell_t> > m_Cells;
  QVector<cell_t*>        m_GpuCells;
  bool                    m_UpdateRequired;


protected: // methods

  void update();


public: // methods

  GPU_CartesianLevelSetBC(PatchGrid* patch_grid, int cuda_device = 0, size_t thread_limit = 0);

  CUDA_DO static void grad(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k, real &gx, real &gy, real &gz);
  CUDA_DO static void getOutsideState(GPU_CartesianPatch &patch,
                                      size_t i, size_t j, size_t k,
                                      int di, int dj, int dk,
                                      real &h0, real &h1, real &h2, real &w, real* var1, real* var2);

  CUDA_HO virtual void operator()();

};


template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::GPU_CartesianLevelSetBC(PatchGrid* patch_grid, int cuda_device, size_t thread_limit)
  : GPU_LevelSetBC<DIM, CartesianPatch, GPU_CartesianPatch>(patch_grid, cuda_device, thread_limit)
{
  m_UpdateRequired = true;
}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::update()
{
  if (!m_UpdateRequired) {
    return;
  }
  m_Cells.resize(this->m_Patches.size());
  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    m_Cells[i_patch].clear();
    int imax = this->m_Patches[i_patch]->sizeI();
    int jmax = this->m_Patches[i_patch]->sizeJ();
    int kmax = this->m_Patches[i_patch]->sizeK();
    for (size_t i = 0; i < imax; ++i) {
      for (size_t j = 0; j < jmax; ++j) {
        for (size_t k = 0; k < kmax; ++k) {
          if (LS::G(*this->m_Patches[i_patch], i, j, k) < 0) {
            bool extrapolate = true;
            int di_1 = -1;
            int di_2 =  1;
            int dj_1 = -1;
            int dj_2 =  1;
            int dk_1 = -1;
            int dk_2 =  1;
            if (i < 2)          di_1 = 0;
            if (i >= imax - 2)  di_2 = 0;
            if (j < 2)          dj_1 = 0;
            if (j >= jmax - 2)  dj_2 = 0;
            if (k < 2)          dk_1 = 0;
            if (k >= kmax - 2)  dk_2 = 0;
            for (int di = di_1; di <= di_2; ++di) {
              for (int dj = dj_1; dj <= dj_2; ++dj) {
                for (int dk = dk_1; dk <= dk_2; ++dk) {
                  if (di != 0 || dj != 0 || dk != 0) {
                    if ((i+2*di) < 0 || (j+2*dj) < 0 || (k+2*dk) < 0) {
                      BUG;
                    }
                    if (!this->m_Patches[i_patch]->checkRange(i, j, k)) {
                      BUG;
                    }
                    if (!this->m_Patches[i_patch]->checkRange(i+2*di, j+2*dj, k+2*dk)) {
                      BUG;
                    }
                    if (LS::G(*this->m_Patches[i_patch], i + di, j + dj, k + dk) >= 0) {
                      extrapolate = false;
                      break;
                    }
                  }
                }
              }
            }
            cell_t cell;
            cell.i = i;
            cell.j = j;
            cell.k = k;
            cell.extrapolate = extrapolate;
            m_Cells[i_patch] << cell;
          }
        }
      }
    }
  }

  // delete old GPU arrays
  foreach (cell_t* cells, m_GpuCells) {
    cudaFree(cells);
    CUDA_CHECK_ERROR;
  }
  m_GpuCells.clear();

  // allocate new arrays
  m_GpuCells.resize(this->m_Patches.size());
  for (int i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    if (m_Cells[i_patch].size() > 0) {
      cudaMalloc(&m_GpuCells[i_patch], m_Cells[i_patch].size()*sizeof(cell_t));
      CUDA_CHECK_ERROR;
      cell_t* cells = new cell_t[m_Cells[i_patch].size()];
      for (int i = 0; i < m_Cells[i_patch].size(); ++i) {
        cells[i] = m_Cells[i_patch][i];
      }
      cudaMemcpy(m_GpuCells[i_patch], cells, m_Cells[i_patch].size()*sizeof(cell_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;
      delete [] cells;
    } else {
      m_GpuCells[i_patch] = NULL;
    }
  }

  m_UpdateRequired = false;
}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::grad(GPU_CartesianPatch &patch, size_t i, size_t j, size_t k, real &gx, real &gy, real &gz)
{
  gx = 0.5*patch.idx()*(LS::G(patch, i+1, j, k) - LS::G(patch, i-1, j, k));
  gy = 0.5*patch.idy()*(LS::G(patch, i, j+1, k) - LS::G(patch, i, j-1, k));
  gz = 0.5*patch.idz()*(LS::G(patch, i, j, k+1) - LS::G(patch, i, j, k-1));
}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::getOutsideState(GPU_CartesianPatch &patch,
                                                                size_t i, size_t j, size_t k,
                                                                int di, int dj, int dk,
                                                                real &h0, real &h1, real &h2, real &w, real *var1, real *var2)
{
  
  // careful with parallel level sets (e.g. flat plate)
  // there might be a 0/0 occurring
  
  dim_t<DIM> dim;

  h0 = LS::G(patch, i, j, k);
  h1 = LS::G(patch, i + di, j + dj, k + dk);
  h2 = LS::G(patch, i + 2*di, j + 2*dj, k + 2*dk);

  patch.getVar(dim, 0, i + di, j + dj, k + dk, var1);
  patch.getVar(dim, 0, i + 2*di, j + 2*dj, k + 2*dk, var2);
  if (h2 < 0) {
    w = 1;
  } else {
    w = h1/(h1 - h0);
  }
}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
__global__ void GPU_CartesianLevelSetBC_kernel(GPU_CartesianPatch patch, typename GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::cell_t* cells, size_t num_cells)
{
  int idx = blockDim.x*blockIdx.x + threadIdx.x;

  if (idx >= num_cells) {
    return;
  }
  size_t i = cells[idx].i;
  size_t j = cells[idx].j;
  size_t k = cells[idx].k;
  bool extrapolate = cells[idx].extrapolate;

  dim_t<DIM> dim;

  if (LS::G(patch, i, j, k) < 0) {

    real var_1[DIM], var_2[DIM];
    real var[DIM], var_bc[DIM];
    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      var[i_var] = 0.0;
    }

    real gx, gy, gz;
    int count = 0;
    real total_weight = 0;
    real h0, h1, h2, w;
    int di_1 = -1;
    int di_2 =  1;
    int dj_1 = -1;
    int dj_2 =  1;
    int dk_1 = -1;
    int dk_2 =  1;
    if (i < 2)                  di_1 = 0;
    if (i >= patch.sizeI() - 2) di_2 = 0;
    if (j < 2)                  dj_1 = 0;
    if (j >= patch.sizeJ() - 2) dj_2 = 0;
    if (k < 2)                  dk_1 = 0;
    if (k >= patch.sizeK() - 2) dk_2 = 0;
    GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::grad(patch, i, j, k, gx, gy, gz);
    for (int di = di_1; di <= di_2; ++di) {
      for (int dj = dj_1; dj <= dj_2; ++dj) {
        for (int dk = dk_1; dk <= dk_2; ++dk) {
          if (di != 0 || dj != 0 || dk != 0) {
            if (LS::G(patch, i + di, j + dj, k + dk) >= 0 || (extrapolate && LS::G(patch, i + di, j + dj, k + dk) > LS::G(patch, i, j, k))) {
              real dx = di*patch.dx();
              real dy = dj*patch.dy();
              real dz = dk*patch.dz();
              real weight = fabs(dx*gx + dy*gy + dz*gz)/sqrt(dx*dx + dy*dy + dz*dz);
              if (extrapolate) {
                patch.getVar(dim, 0, i + di, j + dj, k + dk, var_bc);
              } else {
                GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::getOutsideState(patch, i, j, k, di, dj, dk, h0, h1, h2, w, var_1, var_2);
                BC::operate(var_1, var_2, var_bc, h0, h1, h2, w, gx, gy, gz);
              }
              ++count;
              total_weight += weight;
              for (size_t i_var = 0; i_var < DIM; ++i_var) {
                var[i_var] += weight*var_bc[i_var];
              }
            }
          }
        }
      }
    }
    if (count > 0) {
      patch.getVar(dim, 0, i, j, k, var_bc);
      for (size_t i_var = 0; i_var < DIM - NUM_LS; ++i_var) {
        var_bc[i_var] = var[i_var]/total_weight;
      }
      //real G_self = LS::G(patch, i, j, k);
      patch.setVar(dim, 2, i, j, k, var_bc);
      //LS::updateG(patch, i, j, k, G_self, 2);
    }

  }
}

template <unsigned int DIM, unsigned int NUM_LS, typename LS, typename BC>
void GPU_CartesianLevelSetBC<DIM,NUM_LS,LS,BC>::operator()()
{
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  CUDA_CHECK_ERROR;

  update();
  size_t max_num_threads = this->m_MaxNumThreads;
  cudaDeviceSynchronize();

  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {

    CUDA_CHECK_ERROR;

    if (m_Cells[i_patch].size() > 0) {
      int num_cells   = m_Cells[i_patch].size();
      int num_blocks  = max(int(16), int(num_cells/max_num_threads) + 1);
      int num_threads = num_cells/num_blocks + 1;
      if (num_cells > num_blocks*num_threads) BUG;
      cudaMemcpy(this->m_GpuPatches[i_patch].getField(2), this->m_GpuPatches[i_patch].getField(0), this->m_GpuPatches[i_patch].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
      GPU_CartesianLevelSetBC_kernel<DIM,NUM_LS,LS,BC> <<<num_blocks, num_threads>>>(this->m_GpuPatches[i_patch], m_GpuCells[i_patch], num_cells);
      CUDA_CHECK_ERROR;
      cudaDeviceSynchronize();
      cudaMemcpy(this->m_GpuPatches[i_patch].getField(0), this->m_GpuPatches[i_patch].getField(2), this->m_GpuPatches[i_patch].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
    }

  }
}

#endif // GPU_CARTESIANLEVELSETBC_H
