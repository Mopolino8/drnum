#ifndef GPU_PATCHITERATOR_H
#define GPU_PATCHITERATOR_H

#include "patchiterator.h"
#include "cartesianpatch.h"
#include "gpu_patch.h"

#include <cuda.h>
#include <cuda_runtime_api.h>

class GPU_PatchIterator : public PatchIterator
{

protected:

  GPU_Patch* m_GpuPatch;


public:

  GPU_PatchIterator(Patch* patch) : PatchIterator(patch)
  {
    m_GpuPatch = NULL;
  }

  void copyField(size_t i_src, size_t i_dst)
  {
    if (m_GpuPatch == NULL) {
      BUG;
    }
    cudaMemcpy(m_GpuPatch->getField(i_dst), m_GpuPatch->getField(i_src), m_GpuPatch->fieldSize()*sizeof(real) ,cudaMemcpyDeviceToDevice);
  }

  void gpuInfo()
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
    for (int i = 0; i < count; ++i) {
      cout << endl;
      cout << "device #" << i << ":" << endl;
      cout << "================================================" << endl;
      cudaDeviceProp prop;
      if (cudaGetDeviceProperties(&prop, i) != cudaSuccess) {
        cerr << "error fetching device properties" << endl;
        exit(EXIT_FAILURE);
      }
      cout << "name                   : " << prop.name << endl;
      cout << "warpSize               : " << prop.warpSize << endl;
      cout << "maxThreadsPerBlock     : " << prop.maxThreadsPerBlock << endl;
      cout << "maxThreadsDim[0]       : " << prop.maxThreadsDim[0] << endl;
      cout << "maxThreadsDim[1]       : " << prop.maxThreadsDim[1] << endl;
      cout << "maxThreadsDim[2]       : " << prop.maxThreadsDim[2] << endl;
      cout << "maxGridSize[0]         : " << prop.maxGridSize[0] << endl;
      cout << "maxGridSize[1]         : " << prop.maxGridSize[1] << endl;
      cout << "maxGridSize[2]         : " << prop.maxGridSize[2] << endl;
      cout << "compute capability     : " << prop.major << "." << prop.minor << endl;
      cout << "sharedMemPerBlock [kb] : " << double(prop.sharedMemPerBlock)/1024 << endl;
      cout << "totalGlobalMem [Mb]    : " << double(prop.totalGlobalMem)/(1024*1024) << endl;
      cout << endl;
    }
  }

};

#endif // GPU_PATCHITERATOR_H
