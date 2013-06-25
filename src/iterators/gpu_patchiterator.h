#ifndef GPU_PATCHITERATOR_H
#define GPU_PATCHITERATOR_H

#include "tpatchiterator.h"
#include "gpu_patch.h"
#include "cudatools.h"


template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
class GPU_PatchIterator : public TPatchIterator<T_CPU, DIM, OP>
{

protected: // attributes

  bool          m_GpuPointersSet;
  vector<T_GPU> m_GpuPatches;


public:

  GPU_PatchIterator(OP op);

  CUDA_HO void updateHost();
  CUDA_HO void updateDevice();

  virtual void addPatch(Patch *patch);
  virtual void copyField(size_t i_src, size_t i_dst);
  virtual void copyDonorData(size_t i_field);

};


template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::GPU_PatchIterator(OP op)
  : TPatchIterator<T_CPU, DIM, OP>(op)
{
  m_GpuPointersSet = false;
  //m_GpuPatches.reserve(max(size_t(100), patch_grid.getNumPatches()));
}

template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
void GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::addPatch(Patch *patch)
{
  T_CPU* cpu_patch = dynamic_cast<T_CPU*>(patch);
  if (cpu_patch == NULL) {
    BUG;
  }
  TPatchIterator<T_CPU, DIM, OP>::addPatch(cpu_patch);
  T_GPU gpu_patch(cpu_patch);
  m_GpuPatches.push_back(gpu_patch);
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
  if (!m_GpuPointersSet) {
    for (size_t i = 0; i < m_GpuPatches.size(); ++i) {
      m_GpuPatches[i].updateDonorPointers();
    }
    m_GpuPointersSet = true;
  }
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    m_GpuPatches[i].copyToDevice(this->m_Patches[i]);
  }
}

template <typename T_GPU>
__global__ void GPU_PatchIterator_kernelCopyDonorData(T_GPU patch, size_t i_field)
{
}

template <typename T_CPU, typename T_GPU, unsigned int DIM, typename OP>
void GPU_PatchIterator<T_CPU, T_GPU, DIM, OP>::copyDonorData(size_t i_field)
{
  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    //GPU_PatchIterator_kernelCopyDonorData<<<1,1>>>(this->m_GpuPatches[i_patch], i_field);
  }
}




#endif // GPU_PATCHITERATOR_H
