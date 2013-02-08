#ifndef GPU_CARTESIANPATCH_H
#define GPU_CARTESIANPATCH_H

#include "cartesianpatch.h"
#include <cuda.h>
#include <cuda_runtime_api.h>

#include "gpu_patch.h"

class GPU_CartesianPatch : public GPU_Patch
{

#include "cartesianpatch_common.h"

public:

  CUDA_HO GPU_CartesianPatch(const CartesianPatch& patch)
  {
    copyAttributes(patch);
  }

};



#endif // GPU_CARTESIANPATCH_H
