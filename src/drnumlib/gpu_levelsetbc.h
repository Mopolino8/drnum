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

#ifndef GPU_LEVELSETBC_H
#define GPU_LEVELSETBC_H

template <unsigned int DIM, typename TCPU, typename TGPU>
class GPU_LevelSetBC;

#include "cudatools.h"

#include "drnum.h"
#include "genericoperation.h"

template <unsigned int DIM, typename TCPU, typename TGPU>
class GPU_LevelSetBC : public GenericOperation
{

protected: // attributes

  PatchGrid*     m_PatchGrid;
  vector<TCPU*>  m_Patches;
  vector<TGPU>   m_GpuPatches;
  size_t         m_MaxNumThreads;
  int            m_CudaDevice;


protected: // methods

  CUDA_HO void copyField(size_t i_src, size_t i_dst);


public: // methods

  GPU_LevelSetBC(PatchGrid* patch_grid, int cuda_device = 0, size_t thread_limit = 0);

};


template <unsigned int DIM, typename TCPU, typename TGPU>
GPU_LevelSetBC<DIM,TCPU,TGPU>::GPU_LevelSetBC(PatchGrid* patch_grid, int cuda_device, size_t thread_limit)
{
  m_CudaDevice = cuda_device;
  int count;
  if (cudaGetDeviceCount(&count) != cudaSuccess) {
    cerr << "error detecting CUDA devices" << endl;
    exit(EXIT_FAILURE);
  }
  if (count < m_CudaDevice + 1) {
    CudaTools::info();
    cerr << "specified CUDA device does not exists" << endl;
    exit(EXIT_FAILURE);
  }
  cudaDeviceProp prop;
  if (cudaGetDeviceProperties(&prop, m_CudaDevice) != cudaSuccess) {
    cerr << "error fetching device properties" << endl;
    exit(EXIT_FAILURE);
  }
  cudaSetDevice(m_CudaDevice);
  m_MaxNumThreads = min(prop.maxThreadsPerBlock, prop.maxThreadsDim[0]);
  if (thread_limit > 0) {
    m_MaxNumThreads = min(thread_limit, m_MaxNumThreads);
  }
  m_PatchGrid = patch_grid;
  for (size_t i = 0; i < m_PatchGrid->getNumPatches(); ++i) {
    TCPU* patch = dynamic_cast<TCPU*>(m_PatchGrid->getPatch(i));
    if (patch) {
      m_Patches.push_back(patch);
      m_GpuPatches.push_back(TGPU(patch));
    }
  }
}

template <unsigned int DIM, typename TCPU, typename TGPU>
void GPU_LevelSetBC<DIM, TCPU, TGPU>::copyField(size_t i_src, size_t i_dst)
{
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    cudaMemcpy(m_GpuPatches[i].getField(i_dst), m_GpuPatches[i].getField(i_src), m_GpuPatches[i].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
  }
}

#endif // GPU_LEVELSETBC_H
