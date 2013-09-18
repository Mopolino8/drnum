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
#ifndef GPU_CARTESIANITERATOR_H
#define GPU_CARTESIANITERATOR_H

#include "patchiterator.h"
#include "cartesianpatch.h"
#include "gpu_cartesianpatch.h"
#include "gpu_patchiterator.h"

template <unsigned int DIM, typename OP>
class GPU_CartesianIterator : public GPU_PatchIterator<CartesianPatch, GPU_CartesianPatch<DIM>, OP>
{

public:

  using GPU_PatchIterator<CartesianPatch, GPU_CartesianPatch<DIM>, OP>::addPatch;

  CUDA_HO GPU_CartesianIterator(OP op);

  CUDA_HO virtual void compute(real factor, const vector<size_t>& patches);

};

template <unsigned int DIM, typename OP>
GPU_CartesianIterator<DIM,OP>::GPU_CartesianIterator(OP op)
  : GPU_PatchIterator<CartesianPatch, GPU_CartesianPatch<DIM>, OP>(op)
{
}

#define NDEBUG

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelXFieldFluxes(GPU_CartesianPatch<DIM> patch, OP op, size_t offset)
{
  size_t i = 2*blockIdx.x + offset;
  size_t j = blockDim.y*blockIdx.y + threadIdx.y;
  size_t k = threadIdx.x;
  if (i >= patch.sizeI() || j >= patch.sizeJ()) return;

  real A = patch.dy()*patch.dz();
  real x = (real)0.5*patch.dx() + i*patch.dx();
  real y = (real)0.5*patch.dy() + j*patch.dy();
  real z = (real)0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.xField(&patch, i, j, k, x, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i-1, j, k) -= flux[i_var];
    patch.f(2, i_var, i,   j, k) += flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelYFieldFluxes(GPU_CartesianPatch<DIM> patch, OP op, size_t offset)
{
  size_t i = blockDim.y*blockIdx.x + threadIdx.y;
  size_t j = 2*blockIdx.y + offset;
  size_t k = threadIdx.x;
  if (j >= patch.sizeJ() || i >= patch.sizeI()) return;

  real A = patch.dx()*patch.dz();
  real x = (real)0.5*patch.dx() + i*patch.dx();
  real y = (real)0.5*patch.dy() + j*patch.dy();
  real z = (real)0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.yField(&patch, i, j, k, x, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, j-1, k) -= flux[i_var];
    patch.f(2, i_var, i, j,   k) += flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelZFieldFluxes(GPU_CartesianPatch<DIM> patch, OP op, size_t offset)
{
  size_t i = blockIdx.x;
  size_t j = blockDim.y*blockIdx.y + threadIdx.y;
  size_t k = 2*threadIdx.x + offset;
  if (k >= patch.sizeK() || j >= patch.sizeJ()) return;

  real A = patch.dx()*patch.dy();
  real x = (real)0.5*patch.dx() + i*patch.dx();
  real y = (real)0.5*patch.dy() + j*patch.dy();
  real z = (real)0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.zField(&patch, i, j, k, x, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, j, k-1) -= flux[i_var];
    patch.f(2, i_var, i, j, k  ) += flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelXBoundaryFluxes(GPU_CartesianPatch<DIM> patch, OP op)
{
  size_t j = blockIdx.x;
  size_t k = threadIdx.x;

  real A  = patch.dy()*patch.dz();
  real x1 = (real)0.5*patch.dx();
  real x2 = (real)0.5*patch.dx() + patch.sizeI()*patch.dx();
  real y  = (real)0.5*patch.dy() + j*patch.dy();
  real z  = (real)0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.xWallM(&patch, 0, j, k, x1, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, 0, j, k) += flux[i_var];
  }
  fill(flux, DIM, 0);
  op.xWallP(&patch, patch.sizeI(), j, k, x2, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, patch.sizeI() - 1, j, k) -= flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelYBoundaryFluxes(GPU_CartesianPatch<DIM> patch, OP op)
{
  size_t i = blockIdx.x;
  size_t k = threadIdx.x;

  real A = patch.dx()*patch.dz();
  real x  = (real)0.5*patch.dx() + i*patch.dx();
  real y1 = (real)0.5*patch.dy();
  real y2 = (real)0.5*patch.dy() + patch.sizeJ()*patch.dy();
  real z  = (real)0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.yWallM(&patch, i, 0, k, x, y1, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, 0, k) += flux[i_var];
  }
  fill(flux, DIM, 0);
  op.yWallP(&patch, i, patch.sizeJ(), k, x, y2, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, patch.sizeJ() - 1, k) -= flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelZBoundaryFluxes(GPU_CartesianPatch<DIM> patch, OP op)
{
  size_t i = blockIdx.x;
  size_t j = threadIdx.x;

  real A = patch.dx()*patch.dy();
  real x  = (real)0.5*patch.dx() + i*patch.dx();
  real y  = (real)0.5*patch.dy() + j*patch.dy();
  real z1 = (real)0.5*patch.dz();
  real z2 = (real)0.5*patch.dz() + patch.sizeK()*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.zWallM(&patch, i, j, 0, x, y, z1, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, j, 0) += flux[i_var];
  }
  fill(flux, DIM, 0);
  op.zWallP(&patch, i, j, patch.sizeK(), x, y, z2, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, j, patch.sizeK() - 1) -= flux[i_var];
  }
}

template <unsigned int DIM>
__global__ void GPU_CartesianIterator_kernelAdvance(GPU_CartesianPatch<DIM> patch, real factor)
{
  size_t i = blockIdx.x;
  size_t j = blockIdx.y;
  size_t k = threadIdx.x;

  if (!patch.checkRange(i,j,k)) {
    return;
  }

  factor /= patch.dV();
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(0, i_var, i, j, k) = patch.f(1, i_var, i, j, k) + factor*patch.f(2, i_var, i, j, k);
  }
}

#define LAUNCH_KERNEL(KERNEL, BLOCKS, THREADS, PATCH, OP, OFFSET) \
  KERNEL <<<BLOCKS,THREADS>>> (PATCH, OP, OFFSET);                \
  cudaThreadSynchronize();                                        \
  cudaError_t err = cudaGetLastError();                           \
  if (err != cudaSuccess) {                                       \
    cerr << "\n" << cudaGetErrorString(err) << "\n" << endl;      \
    BUG;                                                          \
  }

template <unsigned int DIM, typename OP>
void GPU_CartesianIterator<DIM,OP>::compute(real factor, const vector<size_t> &patches)
{
  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  CUDA_CHECK_ERROR;

  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {

    cudaMemset(this->m_GpuPatches[i_patch].getField(2), 0, this->m_GpuPatches[i_patch].fieldSize()*sizeof(real));
    CUDA_CHECK_ERROR;

    size_t max_num_threads = GPU_PatchIterator<CartesianPatch, GPU_CartesianPatch<DIM>, OP>::m_MaxNumThreads;
    size_t k_lines = max(size_t(1), size_t(max_num_threads/this->m_Patches[i_patch]->sizeK()));

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI()/2+1, this->m_Patches[i_patch]->sizeJ()/k_lines+1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), k_lines, 1);
      GPU_CartesianIterator_kernelXFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 1);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
      GPU_CartesianIterator_kernelXFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 2);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI()/k_lines+1, this->m_Patches[i_patch]->sizeJ()/2+1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), k_lines, 1);
      GPU_CartesianIterator_kernelYFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 1);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
      GPU_CartesianIterator_kernelYFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 2);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), this->m_Patches[i_patch]->sizeJ()/k_lines+1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK()/2+1, k_lines, 1);
      GPU_CartesianIterator_kernelZFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 1);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
      GPU_CartesianIterator_kernelZFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 2);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
    }


    {
      dim3 blocks(this->m_Patches[i_patch]->sizeJ(), 1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), 1, 1);
      GPU_CartesianIterator_kernelXBoundaryFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), 1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), 1, 1);
      GPU_CartesianIterator_kernelYBoundaryFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), 1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeJ(), 1, 1);
      GPU_CartesianIterator_kernelZBoundaryFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
    }


    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), this->m_Patches[i_patch]->sizeJ(), 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), 1, 1);
      GPU_CartesianIterator_kernelAdvance<DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], factor);
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
    }

  }
}


#endif // GPU_CARTESIANITERATOR_H
