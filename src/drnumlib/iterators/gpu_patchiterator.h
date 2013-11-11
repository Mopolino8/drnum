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
#ifndef GPU_PATCHITERATOR_H
#define GPU_PATCHITERATOR_H

#include "tpatchiterator.h"
#include "gpu_patch.h"
#include "cudatools.h"

#include <set>


template <unsigned int DIM, typename T_CPU, typename T_GPU, typename OP>
class GPU_PatchIterator : public TPatchIterator<T_CPU, OP>
{

protected: // attributes

  bool          m_GpuPointersSet;
  vector<T_GPU> m_GpuPatches;
  size_t        m_MaxNumThreads;
  int           m_CudaDevice;


public:

  GPU_PatchIterator(OP op, int cuda_device = 0, size_t thread_limit = 0);

  CUDA_HO virtual void updateHost();
  CUDA_HO virtual void updateDevice();

  virtual void addPatch(Patch *patch);
  virtual void copyField(size_t i_src, size_t i_dst);
  virtual void copyDonorData(size_t i_field);

};


template <unsigned int DIM, typename T_CPU, typename T_GPU, typename OP>
GPU_PatchIterator<DIM, T_CPU, T_GPU, OP>::GPU_PatchIterator(OP op, int cuda_device, size_t thread_limit)
  : TPatchIterator<T_CPU, OP>(op)
{
  m_GpuPointersSet = false;
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
}

template <unsigned int DIM, typename T_CPU, typename T_GPU, typename OP>
void GPU_PatchIterator<DIM, T_CPU, T_GPU, OP>::addPatch(Patch *patch)
{
  T_CPU* cpu_patch = dynamic_cast<T_CPU*>(patch);
  if (cpu_patch == NULL) {
    BUG;
  }
  TPatchIterator<T_CPU, OP>::addPatch(cpu_patch);
  T_GPU gpu_patch(cpu_patch);
  m_GpuPatches.push_back(gpu_patch);
}

template <unsigned int DIM, typename T_CPU, typename T_GPU, typename OP>
void GPU_PatchIterator<DIM, T_CPU, T_GPU, OP>::copyField(size_t i_src, size_t i_dst)
{
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    cudaMemcpy(m_GpuPatches[i].getField(i_dst), m_GpuPatches[i].getField(i_src), m_GpuPatches[i].fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
  }
}

template <unsigned int DIM, typename T_CPU, typename T_GPU, typename OP>
void GPU_PatchIterator<DIM, T_CPU, T_GPU, OP>::updateHost()
{
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    m_GpuPatches[i].copyFromDevice(this->m_Patches[i]);
  }
}

template <unsigned int DIM, typename T_CPU, typename T_GPU, typename OP>
void GPU_PatchIterator<DIM, T_CPU, T_GPU, OP>::updateDevice()
{
  cout << "void GPU_PatchIterator<T_CPU, T_GPU, OP>::updateDevice()" << endl;
  if (!m_GpuPointersSet) {
    for (size_t i = 0; i < m_GpuPatches.size(); ++i) {
      m_GpuPatches[i].updateDonorPointers(PatchIterator::getPatch(i));
    }
    m_GpuPointersSet = true;
  }
  for (size_t i = 0; i < this->m_Patches.size(); ++i) {
    m_GpuPatches[i].copyToDevice(this->m_Patches[i]);
  }
}

template <unsigned int DIM, typename T_GPU>
__global__ void GPU_PatchIterator_kernelResetReceiverData(T_GPU patch, size_t i_field)
{
  size_t i = blockDim.x*blockIdx.x + threadIdx.x;
  if (i < patch.getNumReceivingCellsUnique()) {
    size_t i_rec = patch.getReceivingCellIndicesUnique()[i];
    for (size_t i_var = 0; i_var < DIM; ++i_var) {
      patch.getVariable(i_field, i_var)[i_rec] = 0;
    }
  }
}

template <unsigned int DIM, unsigned int STRIDE, typename T_GPU>
__global__ void GPU_PatchIterator_kernelCopyDonorData(T_GPU patch, size_t i_field, size_t i_donor)
{
  size_t i = blockDim.x*blockIdx.x + threadIdx.x;
  donor_t donor = patch.getDonors()[i_donor];
  if (i < donor.num_receiver_cells) {

    // receiving cells index
    size_t i_rec = patch.getReceivingCellIndicesConcat()[donor.receiver_index_field_start + i];

    // start address in m_DonorCells/m_DonorWeights pattern
    size_t i_donor_cells_start = donor.donor_wi_field_start + i*donor.stride;

    // loop for contributing cells
    for (size_t i_contrib = 0; i_contrib < STRIDE; ++i_contrib) {

      size_t i_wi              = i_donor_cells_start + i_contrib;      // index of donor cell in concatenated lists
      size_t donor_cell_index  = patch.getDonorIndexConcat()[i_wi];
      real   donor_cell_weight = patch.getDonorWeightConcat()[i_wi];

      for (size_t i_var = 0; i_var < DIM; ++i_var) {
        real* dvar = donor.data + i_var*donor.variable_size;
        //donated_var[i_var] = dvar[donor_cell_index];
        patch.getVariable(i_field, i_var)[i_rec] += donor_cell_weight*dvar[donor_cell_index];
        //patch.getVariable(i_field, i_var)[i_rec] += 1.0;//donor_cell_weight;
      }

      /*
      // transform vector variables
      for (size_t i_vec = 0; i_vec < patch.getNumVectorVars(); ++i_vec) {
        size_t i_var = patch.getVectorVarIndices()[i_vec];
        real u = donor.axx*donated_var[i_var + 0] + donor.axy*donated_var[i_var + 1] + donor.axz*donated_var[i_var + 2];
        real v = donor.ayx*donated_var[i_var + 0] + donor.ayy*donated_var[i_var + 1] + donor.ayz*donated_var[i_var + 2];
        real w = donor.azx*donated_var[i_var + 0] + donor.azy*donated_var[i_var + 1] + donor.azz*donated_var[i_var + 2];
        donated_var[i_var + 0] = u;
        donated_var[i_var + 1] = v;
        donated_var[i_var + 2] = w;
      }

      // contribute to receiving cell
      for (size_t i_var = 0; i_var < patch.numVariables(); ++i_var) {
        patch.getVariable(i_field, i_var)[i_rec] += donor_cell_weight*donated_var[i_var];
      }

      delete [] donated_var;
      */
    }
  }
}

template <unsigned int DIM, typename T_CPU, typename T_GPU, typename OP>
void GPU_PatchIterator<DIM, T_CPU, T_GPU, OP>::copyDonorData(size_t i_field)
{
  // reset receiver cells
  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    int N = PatchIterator::getPatch(i_patch)->getNumReceivingCellsUnique();
    int blocks  = N/m_MaxNumThreads + 1;
    GPU_PatchIterator_kernelResetReceiverData <DIM, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>>(m_GpuPatches[i_patch], i_field);
    CUDA_CHECK_ERROR
  }

  cudaThreadSynchronize();
  CUDA_CHECK_ERROR;

  // compute interpolated data
  for (size_t i_patch = 0; i_patch < this->m_Patches.size(); ++i_patch) {
    for (size_t i_donor = 0; i_donor < m_GpuPatches[i_patch].getNumDonorPatches(); ++i_donor) {
      int N = PatchIterator::getPatch(i_patch)->getDonors()[i_donor].num_receiver_cells;
      int stride = PatchIterator::getPatch(i_patch)->getDonors()[i_donor].stride;
      int blocks  = N/m_MaxNumThreads + 1;
      if        (stride == 1) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 1, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else if (stride == 2) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 2, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else if (stride == 3) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 3, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else if (stride == 4) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 4, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else if (stride == 5) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 5, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else if (stride == 6) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 6, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else if (stride == 7) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 7, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else if (stride == 8) {
        GPU_PatchIterator_kernelCopyDonorData <DIM, 8, GPU_CartesianPatch> <<<blocks, m_MaxNumThreads>>> (m_GpuPatches[i_patch], i_field, i_donor);
      } else {
        BUG;
      }
      CUDA_CHECK_ERROR;
      //cudaThreadSynchronize();
    }
  }
}




#endif // GPU_PATCHITERATOR_H
