#ifndef GPU_PATCH_H
#define GPU_PATCH_H

#include "patch.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

class GPU_Patch
{

#include "patch_common.h"

public:

  CUDA_HO GPU_Patch(const Patch& patch)
  {
    copyAttributes(patch);
    cudaMalloc(&m_Data, sizeof(real)*patch.dataSize());
    copyToDevice(patch);
  }

  CUDA_HO copyToDevice(const Patch& patch)
  {
    cudaMemcpy(m_Data, patch.getData(), patch.dataSize(), cudaMemcpyHostToDevice);
  }

  CUDA_HO copyFromDevice(Patch& patch)
  {
    cudaMemcpy(patch.getData(), m_Data, dataSize(), cudaMemcpyDeviceToHost);
  }

};

#endif // GPU_PATCH_H
