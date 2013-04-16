#ifndef GPU_PATCHITERATOR_H
#define GPU_PATCHITERATOR_H

#include "tpatchiterator.h"
#include "gpu_patch.h"
#include "cudatools.h"

template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
class GPU_PatchIterator : public TPatchIterator<T_CPU, DIM, OP>
{

protected: // attributes

  vector<T_GPU> m_GpuPatches;


public:

  GPU_PatchIterator(PatchGrid &patch_grid, OP op);

  CUDA_HO void addPatch(T_CPU* patch);
  CUDA_HO void updateHost();
  CUDA_HO void updateDevice();

  virtual void copyField(size_t i_src, size_t i_dst);

};


template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::GPU_PatchIterator(PatchGrid &patch_grid, OP op)
  : TPatchIterator<T_CPU, DIM, OP>(patch_grid, op)
{
  m_GpuPatches.reserve(max(size_t(100), patch_grid.getNumPatches()));
}

template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
void GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::addPatch(T_CPU *patch)
{
  TPatchIterator<T_CPU, DIM, OP>::addPatch(patch);
  m_GpuPatches.push_back(T_GPU(patch));
}

template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
void GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::copyField(size_t i_src, size_t i_dst)
{
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    cudaMemcpy(m_GpuPatches[i].getField(i_dst), m_GpuPatches[i].getField(i_src), m_GpuPatches[i].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
  }
}

template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
void GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::updateHost()
{
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    m_GpuPatches[i].copyFromDevice(this->m_Patches[i]);
  }
}

template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
void GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::updateDevice()
{
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    m_GpuPatches[i].copyToDevice(this->m_Patches[i]);
  }
}



#endif // GPU_PATCHITERATOR_H