#ifndef GPU_CARTESIANITERATOR_H
#define GPU_CARTESIANITERATOR_H

#include "patchiterator.h"
#include "cartesianpatch.h"
#include "gpu_cartesianpatch.h"
#include "gpu_patchiterator.h"

template <unsigned int DIM, typename OP>
class GPU_CartesianIterator : public GPU_PatchIterator<CartesianPatch, GPU_CartesianPatch, DIM, OP>
{

private: // attributes

  size_t m_MaxNumThreads;
  //size_t m_NumBlocks;
  //size_t m_NumThreads;


protected: // attributes


public:

  using GPU_PatchIterator<CartesianPatch, GPU_CartesianPatch, DIM, OP>::addPatch;

  CUDA_HO GPU_CartesianIterator(PatchGrid &patch_grid, OP op);

  CUDA_HO virtual void compute(real factor, const vector<size_t>& patches);

};

template <unsigned int DIM, typename OP>
GPU_CartesianIterator<DIM,OP>::GPU_CartesianIterator(PatchGrid &patch_grid, OP op)
  : GPU_PatchIterator<CartesianPatch, GPU_CartesianPatch, DIM, OP>(patch_grid, op)
{
  int count;
  if (cudaGetDeviceCount(&count) != cudaSuccess) {
    cerr << "error detecting CUDA devices" << endl;
    exit(EXIT_FAILURE);
  }
  if (count == 0) {
    cerr << "no CUDA devices found" << endl;
    exit(EXIT_FAILURE);
  }
  cudaDeviceProp prop;
  if (cudaGetDeviceProperties(&prop, 0) != cudaSuccess) {
    cerr << "error fetching device properties" << endl;
    exit(EXIT_FAILURE);
  }
  m_MaxNumThreads = min(10000,min(prop.maxThreadsPerBlock, prop.maxThreadsDim[0]));

  //m_NumThreads = patch_grid->sizeK();//128;
  //m_NumBlocks = patch_grid->variableSize()/m_NumThreads + 1;
}


template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelXFieldFluxes(GPU_CartesianPatch patch, OP op, size_t offset)
{
  size_t i = 2*blockIdx.x + offset;
  size_t j = blockDim.y*blockIdx.y + threadIdx.y;
  size_t k = threadIdx.x;
  if (i >= patch.sizeI() || j >= patch.sizeJ()) return;

  real A = patch.dy()*patch.dz();
  real x = 0.5*patch.dx() + i*patch.dx();
  real y = 0.5*patch.dy() + j*patch.dy();
  real z = 0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.xField(&patch, i, j, k, x, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i-1, j, k) -= flux[i_var];
    patch.f(2, i_var, i,   j, k) += flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelYFieldFluxes(GPU_CartesianPatch patch, OP op, size_t offset)
{
  size_t i = blockDim.y*blockIdx.x + threadIdx.y;
  size_t j = 2*blockIdx.y + offset;
  size_t k = threadIdx.x;
  if (j >= patch.sizeJ() || i >= patch.sizeI()) return;

  real A = patch.dx()*patch.dz();
  real x = 0.5*patch.dx() + i*patch.dx();
  real y = 0.5*patch.dy() + j*patch.dy();
  real z = 0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.yField(&patch, i, j, k, x, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, j-1, k) -= flux[i_var];
    patch.f(2, i_var, i, j,   k) += flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelZFieldFluxes(GPU_CartesianPatch patch, OP op, size_t offset)
{
  size_t i = blockIdx.x;
  size_t j = blockDim.y*blockIdx.y + threadIdx.y;
  size_t k = 2*threadIdx.x + offset;
  if (k >= patch.sizeK() || j >= patch.sizeJ()) return;

  real A = patch.dx()*patch.dy();
  real x = 0.5*patch.dx() + i*patch.dx();
  real y = 0.5*patch.dy() + j*patch.dy();
  real z = 0.5*patch.dz() + k*patch.dz();
  real flux[DIM];
  fill(flux, DIM, 0);
  op.zField(&patch, i, j, k, x, y, z, A, flux);
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    patch.f(2, i_var, i, j, k-1) -= flux[i_var];
    patch.f(2, i_var, i, j, k  ) += flux[i_var];
  }
}

template <unsigned int DIM, typename OP>
__global__ void GPU_CartesianIterator_kernelXBoundaryFluxes(GPU_CartesianPatch patch, OP op)
{
  size_t j = blockIdx.x;
  size_t k = threadIdx.x;

  real A  = patch.dy()*patch.dz();
  real x1 = 0.5*patch.dx();
  real x2 = 0.5*patch.dx() + patch.sizeI()*patch.dx();
  real y  = 0.5*patch.dy() + j*patch.dy();
  real z  = 0.5*patch.dz() + k*patch.dz();
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
__global__ void GPU_CartesianIterator_kernelYBoundaryFluxes(GPU_CartesianPatch patch, OP op)
{
  size_t i = blockIdx.x;
  size_t k = threadIdx.x;

  real A = patch.dx()*patch.dz();
  real x  = 0.5*patch.dx() + i*patch.dx();
  real y1 = 0.5*patch.dy();
  real y2 = 0.5*patch.dy() + patch.sizeJ()*patch.dy();
  real z  = 0.5*patch.dz() + k*patch.dz();
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
__global__ void GPU_CartesianIterator_kernelZBoundaryFluxes(GPU_CartesianPatch patch, OP op)
{
  size_t i = blockIdx.x;
  size_t j = threadIdx.x;

  real A = patch.dx()*patch.dy();
  real x  = 0.5*patch.dx() + i*patch.dx();
  real y  = 0.5*patch.dy() + j*patch.dy();
  real z1 = 0.5*patch.dz();
  real z2 = 0.5*patch.dz() + patch.sizeK()*patch.dz();
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
__global__ void GPU_CartesianIterator_kernelAdvance(GPU_CartesianPatch patch, real factor)
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

template <unsigned int DIM, typename OP>
void GPU_CartesianIterator<DIM,OP>::compute(real factor, const vector<size_t> &patches)
{

#ifdef DEBUG
  BUG;
#endif

  cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
  CudaTools::checkError();

  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {

    cudaMemset(this->m_GpuPatches[i_patch].getField(2), 0, this->m_GpuPatches[i_patch].fieldSize()*sizeof(real));
    CudaTools::checkError();

    size_t k_lines = max(size_t(1), size_t(512/this->m_Patches[i_patch]->sizeK()));

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI()/2+1, this->m_Patches[i_patch]->sizeJ()/k_lines+1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), k_lines, 1);
      GPU_CartesianIterator_kernelXFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 1);
      GPU_CartesianIterator_kernelXFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 2);
      CudaTools::checkError();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI()/k_lines+1, this->m_Patches[i_patch]->sizeJ()/2+1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), k_lines, 1);
      GPU_CartesianIterator_kernelYFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 1);
      GPU_CartesianIterator_kernelYFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 2);
      CudaTools::checkError();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), this->m_Patches[i_patch]->sizeJ()/k_lines+1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK()/2+1, k_lines, 1);
      GPU_CartesianIterator_kernelZFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 1);
      GPU_CartesianIterator_kernelZFieldFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op, 2);
      CudaTools::checkError();
    }


    {
      dim3 blocks(this->m_Patches[i_patch]->sizeJ(), 1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), 1, 1);
      GPU_CartesianIterator_kernelXBoundaryFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op);
      CudaTools::checkError();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), 1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), 1, 1);
      GPU_CartesianIterator_kernelYBoundaryFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op);
      CudaTools::checkError();
    }

    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), 1, 1);
      dim3 threads(this->m_Patches[i_patch]->sizeJ(), 1, 1);
      GPU_CartesianIterator_kernelZBoundaryFluxes <DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], this->m_Op);
      CudaTools::checkError();
    }


    {
      dim3 blocks(this->m_Patches[i_patch]->sizeI(), this->m_Patches[i_patch]->sizeJ(), 1);
      dim3 threads(this->m_Patches[i_patch]->sizeK(), 1, 1);
      GPU_CartesianIterator_kernelAdvance<DIM> <<<blocks, threads>>>(this->m_GpuPatches[i_patch], factor);
      CudaTools::checkError();
    }

  }
}


#endif // GPU_CARTESIANITERATOR_H
