#ifndef GPU_PATCH_H
#define GPU_PATCH_H

#include "patch.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

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
  }

  CUDA_HO void copyToDevice(Patch* patch)
  {
    cudaMemcpy(m_Data, patch->getData(), patch->dataSize()*sizeof(real), cudaMemcpyHostToDevice);
  }

  CUDA_HO void copyFromDevice(Patch* patch)
  {
    cudaMemcpy(patch->getData(), m_Data, dataSize()*sizeof(real), cudaMemcpyDeviceToHost);
  }

};

#endif // GPU_PATCH_H
