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
#ifndef GPU_CARTESIANPATCH_H
#define GPU_CARTESIANPATCH_H

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "cartesianpatch.h"
#include "gpu_patch.h"

class GPU_CartesianPatch;

template <unsigned int DIM>
__global__ void GPU_CartesianPatch_kernelPartialCopy(GPU_CartesianPatch patch, size_t i_field, real *tmp_array,
                                                     size_t i1, size_t i2,
                                                     size_t j1, size_t j2,
                                                     size_t k1, size_t k2);

class GPU_CartesianPatch : public GPU_Patch
{

#include "cartesianpatch_common.h"

public:

  using GPU_Patch::getVar;

  CUDA_HO GPU_CartesianPatch(CartesianPatch* patch) : GPU_Patch(patch)
  {
    copyAttributes(patch);
  }

  template <unsigned int DIM>
  CUDA_HO void partialCopyFromDevice(CartesianPatch *patch, size_t i_field,
                                     size_t i1, size_t i2,
                                     size_t j1, size_t j2,
                                     size_t k1, size_t k2)
  {
    CUDA_CHECK_ERROR
    size_t dim = patch->numVariables();
    size_t tmp_array_length = DIM*(i2-i1)*(j2-j1)*(k2-k1);
    real *tmp_device_array, *tmp_host_array;
    cudaThreadSynchronize();
    cudaMalloc(&tmp_device_array, tmp_array_length*sizeof(real));
    CUDA_CHECK_ERROR
    tmp_host_array = new real[tmp_array_length];
    GPU_CartesianPatch_kernelPartialCopy <DIM> <<<1, 1>>> (*this, i_field, tmp_device_array, i1, i2, j1, j2, k1, k2);
    cudaThreadSynchronize();
    CUDA_CHECK_ERROR
    cudaMemcpy(tmp_host_array, tmp_device_array, tmp_array_length*sizeof(real), cudaMemcpyDeviceToHost);
    CUDA_CHECK_ERROR
    size_t i_tmp = 0;
    for (size_t i_var = 0; i_var < dim; ++i_var) {
      for (size_t i = i1; i < i2; ++i) {
        for (size_t j = j1; j < j2; ++j) {
          for (size_t k = k1; k < k2; ++k) {
            patch->f(i_field, i_var, i, j, k) = tmp_host_array[i_tmp];
            ++i_tmp;
          }
        }
      }
    }
    cudaFree(tmp_device_array);
    CUDA_CHECK_ERROR
    delete [] tmp_host_array;
  }

};


template <unsigned int DIM>
__global__ void GPU_CartesianPatch_kernelPartialCopy(GPU_CartesianPatch patch, size_t i_field, real *tmp_array,
                                                     size_t i1, size_t i2,
                                                     size_t j1, size_t j2,
                                                     size_t k1, size_t k2)
{
  size_t i_tmp = 0;
  for (size_t i_var = 0; i_var < DIM; ++i_var) {
    for (size_t i = i1; i < i2; ++i) {
      for (size_t j = j1; j < j2; ++j) {
        for (size_t k = k1; k < k2; ++k) {
          tmp_array[i_tmp] = patch.f(i_field, i_var, i, j, k);
          ++i_tmp;
        }
      }
    }
  }
}



#endif // GPU_CARTESIANPATCH_H
