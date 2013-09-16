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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
