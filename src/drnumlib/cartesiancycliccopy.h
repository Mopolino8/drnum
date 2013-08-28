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
#ifndef CARTESIANCYCLICCOPY_H
#define CARTESIANCYCLICCOPY_H

class CartesianCyclicCopy;

#include "patchgrid.h"
#include "genericoperation.h"
#include "cartesianpatch.h"
#include "gpu_cartesianpatch.h"

class CartesianCyclicCopy : public GenericOperation
{

private:

  PatchGrid *m_PatchGrid;
  vector<GPU_CartesianPatch*> m_GpuPatches;

  void updatePatches();

public:

  CartesianCyclicCopy(PatchGrid *patch_grid);
  virtual void operator()();

};

#ifdef CUDA
__global__ void CartesianCyclicCopy_kernel(GPU_CartesianPatch patch)
{
  size_t i = blockIdx.x;
  size_t j = threadIdx.x;
  {
    size_t k_from = 1;
    size_t k_to   = patch.sizeK() - 1;
    for (size_t i_var = 0; i_var < patch.numVariables(); ++i_var) {
      patch.f(0, i_var, i, j, k_to) = patch.f(0, i_var, i, j, k_from);
    }
  }
  {
    size_t k_from = patch.sizeK() - 2;
    size_t k_to   = 0;
    for (size_t i_var = 0; i_var < patch.numVariables(); ++i_var) {
      patch.f(0, i_var, i, j, k_to) = patch.f(0, i_var, i, j, k_from);
    }
  }
}
#endif


inline CartesianCyclicCopy::CartesianCyclicCopy(PatchGrid *patch_grid)
{
  m_PatchGrid = patch_grid;
  updatePatches();
}

inline void CartesianCyclicCopy::updatePatches()
{
  if (m_GpuPatches.size() != m_PatchGrid->getNumPatches()) {
    m_GpuPatches.resize(m_PatchGrid->getNumPatches(), NULL);
    for (int i = 0; i < m_PatchGrid->getNumPatches(); ++i) {
      CartesianPatch *cart_patch = dynamic_cast<CartesianPatch*>(m_PatchGrid->getPatch(i));
      if (cart_patch) {
        m_GpuPatches[i] = new GPU_CartesianPatch(cart_patch);
      }
    }
  }
}

inline void CartesianCyclicCopy::operator ()()
{
  updatePatches();
  for (int i = 0; i < m_GpuPatches.size(); ++i) {
    if (m_GpuPatches[i]) {

#ifdef CUDA
      dim3 blocks(this->m_GpuPatches[i]->sizeI(), 1, 1);
      dim3 threads(this->m_GpuPatches[i]->sizeJ(), 1, 1);
      CartesianCyclicCopy_kernel<<<blocks, threads>>>(*(this->m_GpuPatches[i]));
      CUDA_CHECK_ERROR;
      cudaThreadSynchronize();
#else
      BUG;
#endif

    }
  }
}

#endif // CARTESIANCYCLICCOPY_H
