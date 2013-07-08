#ifndef GPU_CARTESIANPATCH_H
#define GPU_CARTESIANPATCH_H

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "cartesianpatch.h"
#include "gpu_patch.h"

class GPU_CartesianPatch : public GPU_Patch
{

#include "cartesianpatch_common.h"

public:

  CUDA_HO GPU_CartesianPatch(CartesianPatch* patch) : GPU_Patch(patch)
  {
    copyAttributes(patch);
  }

};



#endif // GPU_CARTESIANPATCH_H
