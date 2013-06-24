#ifndef GPU_PATCH_H
#define GPU_PATCH_H

#include "patch.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

#include <unordered_map>

class GPU_Patch
{

#include "patch_common.h"

public:

  CUDA_HO GPU_Patch(Patch* patch)
  {
    copyAttributes(patch);
    if (cudaMalloc(&m_Data, sizeof(real)*patch->dataSize()) != cudaSuccess) {
      BUG;
    } else {
      cout << real(sizeof(real)*patch->dataSize())/1024/1024 << " Mbytes allocated on the GPU"<< endl;
    }
    cout << "GPU_Patch::GPU_Patch: m_Data = " << m_Data << endl;
    copyToDevice(patch);

    if (cudaMalloc(&m_ReceivingCellIndicesConcat, sizeof(size_t)*m_NumReceivingCellsConcat) != cudaSuccess) {
      BUG;
    }
    cudaMemcpy(m_ReceivingCellIndicesConcat, patch->getReceivingCellIndicesConcat(), m_NumReceivingCellsConcat*sizeof(size_t), cudaMemcpyHostToDevice);

    if (cudaMalloc(&m_ReceivingCellIndicesUnique, sizeof(size_t)*m_NumReceivingCellsUnique) != cudaSuccess) {
      BUG;
    }
    cudaMemcpy(m_ReceivingCellIndicesUnique, patch->getReceivingCellIndicesUnique(), m_NumReceivingCellsUnique*sizeof(size_t), cudaMemcpyHostToDevice);

    if (cudaMalloc(&m_DonorIndexConcat, sizeof(size_t)*m_NumDonorWIConcat) != cudaSuccess) {
      BUG;
    }
    cudaMemcpy(m_DonorIndexConcat, patch->getDonorIndexConcat(), m_NumDonorWIConcat*sizeof(size_t), cudaMemcpyHostToDevice);

    if (cudaMalloc(&m_DonorWeightConcat, sizeof(real)*m_NumDonorWIConcat) != cudaSuccess) {
      BUG;
    }
    cudaMemcpy(m_DonorWeightConcat, patch->getDonorWeightConcat(), m_NumDonorWIConcat*sizeof(real), cudaMemcpyHostToDevice);

    if (cudaMalloc(&m_Donors, sizeof(donor_t)*m_NumDonorPatches) != cudaSuccess) {
      BUG;
    }
    cudaMemcpy(m_Donors, patch->getDonors(), m_NumDonors*sizeof(donor_t), cudaMemcpyHostToDevice);

  }

  CUDA_HO void copyToDevice(Patch* patch)
  {
    cudaMemcpy(m_Data, patch->getData(), patch->dataSize()*sizeof(real), cudaMemcpyHostToDevice);
  }

  CUDA_HO void copyFromDevice(Patch* patch)
  {
    cudaMemcpy(patch->getData(), m_Data, dataSize()*sizeof(real), cudaMemcpyDeviceToHost);
  }

  CUDA_HO void updateDonorPointers(std::unordered_map<real*,real*> cpu2gpu)
  {
    for (size_t i = 0; i < m_NumDonorPatches; ++i) {
      m_Donors[i].data = (*cpu2gpu)[m_Donors[i].data];
    }
  }

};

#endif // GPU_PATCH_H
