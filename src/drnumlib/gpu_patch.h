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
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef GPU_PATCH_H
#define GPU_PATCH_H

#include "patch.h"
#include "cudatools.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

class GPU_Patch
{

#include "patch_common.h"

protected:

  size_t  m_NumVectorVars;     ///< number of vectorial variables
  size_t* m_VectorVarIndices;  ///< array containing the starting indices of all vectorial variables (e.g. velocity, momentum, ...)


public:

  CUDA_HO GPU_Patch(Patch* patch)
  {
    copyAttributes(patch);
    if (patch->m_GpuDataSet) {
      m_Data                       = patch->m_GpuData;
      m_DonorIndexConcat           = patch->m_GpuDonorIndexConcat;
      m_Donors                     = patch->m_GpuDonors;
      m_DonorWeightConcat          = patch->m_GpuDonorWeightConcat;
      m_ReceivingCellIndicesConcat = patch->m_GpuReceivingCellIndicesConcat;
      m_ReceivingCellIndicesUnique = patch->m_GpuReceivingCellIndicesUnique;
    } else {
      if (cudaMalloc(&m_Data, sizeof(real)*patch->dataSize()) != cudaSuccess) {
        BUG;
      }
      copyToDevice(patch);

      cudaMalloc(&m_ReceivingCellIndicesConcat, sizeof(size_t)*m_NumReceivingCellsConcat);
      CUDA_CHECK_ERROR;
      cudaMemcpy(m_ReceivingCellIndicesConcat, patch->getReceivingCellIndicesConcat(), m_NumReceivingCellsConcat*sizeof(size_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;

      cudaMalloc(&m_ReceivingCellIndicesUnique, sizeof(size_t)*m_NumReceivingCellsUnique);
      CUDA_CHECK_ERROR;
      cudaMemcpy(m_ReceivingCellIndicesUnique, patch->getReceivingCellIndicesUnique(), m_NumReceivingCellsUnique*sizeof(size_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;

      cudaMalloc(&m_DonorIndexConcat, sizeof(size_t)*m_NumDonorWIConcat);
      CUDA_CHECK_ERROR;
      cudaMemcpy(m_DonorIndexConcat, patch->getDonorIndexConcat(), m_NumDonorWIConcat*sizeof(size_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;

      cudaMalloc(&m_DonorWeightConcat, sizeof(real)*m_NumDonorWIConcat);
      CUDA_CHECK_ERROR;
      cudaMemcpy(m_DonorWeightConcat, patch->getDonorWeightConcat(), m_NumDonorWIConcat*sizeof(real), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;

      cudaMalloc(&m_Donors, sizeof(donor_t)*m_NumDonorPatches);
      CUDA_CHECK_ERROR;
      cudaMemcpy(m_Donors, patch->getDonors(), m_NumDonorPatches*sizeof(donor_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;

      // copy data for vectorial variables
      vector<size_t> vector_indices = patch->getVectorVarIndices();
      m_NumVectorVars = vector_indices.size();
      cudaMalloc(&m_VectorVarIndices, sizeof(size_t)*m_NumVectorVars);
      CUDA_CHECK_ERROR;
      cudaMemcpy(m_VectorVarIndices, &vector_indices[0], m_NumVectorVars*sizeof(size_t), cudaMemcpyHostToDevice);
      CUDA_CHECK_ERROR;

      patch->m_GpuData                       = m_Data;
      patch->m_GpuDonorIndexConcat           = m_DonorIndexConcat;
      patch->m_GpuDonors                     = m_Donors;
      patch->m_GpuDonorWeightConcat          = m_DonorWeightConcat;
      patch->m_GpuReceivingCellIndicesConcat = m_ReceivingCellIndicesConcat;
      patch->m_GpuReceivingCellIndicesUnique = m_ReceivingCellIndicesUnique;
      patch->m_GpuDataSet = true;
    }
  }

  CUDA_HO void copyToDevice(Patch* patch)
  {
    cudaMemcpy(m_Data, patch->getData(), patch->dataSize()*sizeof(real), cudaMemcpyHostToDevice);
  }

  CUDA_HO void copyFromDevice(Patch* patch)
  {
    cudaMemcpy(patch->getData(), m_Data, dataSize()*sizeof(real), cudaMemcpyDeviceToHost);
  }

  CUDA_HO void updateDonorPointers(Patch* patch)
  {    
    // get donor array from the GPU
    donor_t* donors = new donor_t [getNumDonorPatches()];
    cudaMemcpy(donors, m_Donors, getNumDonorPatches()*sizeof(donor_t), cudaMemcpyDeviceToHost);

    size_t N = patch->accessNumNeighbours();
    real** cpu = new real* [N];
    real** gpu = new real* [N];
    for (size_t i = 0; i < N; ++i) {
      cpu[i] = patch->accessNeighbour(i)->getData();
      gpu[i] = patch->accessNeighbour(i)->getGpuData();
    }
    for (size_t i = 0; i < m_NumDonorPatches; ++i) {
      bool found = false;
      for (size_t j = 0; j < N; ++j) {
        if (cpu[j] == donors[i].data) {
          found = true;
          donors[i].data = gpu[j];
          break;
        }
      }
      if (!found) {
        BUG;
      }
    }

    // write donor array back to the GPU
    cudaMemcpy(m_Donors, donors, getNumDonorPatches()*sizeof(donor_t), cudaMemcpyHostToDevice);

    delete [] donors;
  }

  CUDA_DH size_t getNumVectorVars()
  {
    return m_NumVectorVars;
  }

  CUDA_DO size_t* getVectorVarIndices()
  {
    return m_VectorVarIndices;
  }

};

#endif // GPU_PATCH_H
