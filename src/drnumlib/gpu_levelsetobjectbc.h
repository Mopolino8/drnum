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
#ifndef GPU_LEVELSETOBJECTBC_H
#define GPU_LEVELSETOBJECTBC_H

template <typename OP>
class GPU_LevelSetObjectBC;

#include "levelsetobjectbc.h"
#include "lslayerdataextrapol.h"

template <typename OP>
class GPU_LevelSetObjectBC : public LevelSetObjectBC
{

private:

  OP m_Op;
  vector<GPU_Patch*> m_GpuPatches;
  LSLayerDataExtrapol* m_TransInnerCellsLayers;
  LSLayerDataExtrapol* m_GpuInnerCellsLayers;
  //  LSLayerDataExtrapol* m_TransOuterCellsLayers; // not needed at present
  //  LSLayerDataExtrapol* m_GpuOuterCellsLayers;

  void update();

public:
  GPU_LevelSetObjectBC(size_t field,
                       LevelSetObject* lso,
                       size_t abuse_field,
                       OP op);

  virtual void operator ()();

};


#ifdef CUDA
template <typename OP>
__global__ void GPU_LevelSetObjectBC_access_kernel(OP op,
                                                   size_t field,
                                                   size_t abuse_field,
                                                   ???m_GpuInnerCellsLayers???) //wie bekommt man ein array hier rein??
{
  ???

  m_Op.access(m_GpuInnerCellsLayers[ll_c], field, abuse_field);

  ???
}


template <typename OP>
__global__ void GPU_LevelSetObjectBC_operateInner_kernel(OP op,
                                                         size_t field,
                                                         size_t abuse_field,
                                                         real relax,
                                                         ???m_GpuInnerCellsLayers???) //wie bekommt man ein array hier rein??
{
  ???

  m_Op.operateInner(Gm_GpuInnerCellsLayers[ll_c], field, abuse_field);

  ???
}


#endif


template <typename OP>
GPU_LevelSetObjectBC<OP>::GPU_LevelSetObjectBC(size_t field,
                                               LevelSetObject* lso,
                                               size_t abuse_field,
                                               OP op)
{
  m_Op = op;
  update();
}


template <typename OP>
void GPU_LevelSetObjectBC<OP>::update()
{
  // get access to gpu patches
  if (m_GpuPatches.size() != m_PatchGrid->getNumPatches()) {
    m_GpuPatches.resize(m_PatchGrid->getNumPatches(), NULL);
    for (int i = 0; i < m_PatchGrid->getNumPatches(); ++i) {
      Patch* patch = m_PatchGrid->getPatch(i);
      m_GpuPatches[i] = new GPU_Patch(patch);
    }
  }
  // build up gpu layer lists
  //.. Transfer Lists: replace m_Data pointer in lists
  m_TransInnerCellsLayers = new LSLayerDataExtrapol[m_InnerCellsLayers->size()];
  //.... Inner
  size_t num_patches = m_PatchGrid->getNumPatches();
  Patch* pivot_patch = m_PatchGrid->getPatch(0);
  for (size_t ll_c = 0; ll_c < m_NumInnerLayerCells; ll_c++) {
    //.... copy all content
    m_TransInnerCellsLayers[ll_c] = m_InnerCellsLayers[ll_c];
    //.... replace data pointer
    if(m_InnerCellsLayers[ll_c].m_Data == pivot_patch->getData()) {
      m_TransInnerCellsLayers[ll_c].m_Data = pivot_patch->getGpuData();
    } else {
      bool found = false;
      for(size_t i_p = 0; i_p < num_patches; i_p++) {
        Patch* patch = m_PatchGrid->getPatch(i_p);
        if(m_InnerCellsLayers[ll_c].m_Data == patch->getData()) {
          m_TransInnerCellsLayers[ll_c].m_Data = patch->getData();
          pivot_patch = patch;
          found = true;
          break;
        }
      }
      if(!found) {BUG;}
    }
  }
  //  //.... Outer
  //  m_TransOuterCellsLayers = new LSLayerDataExtrapol[m_OuterCellsLayers->size()];
  //  pivot_patch = m_PatchGrid->getPatch(0);
  //  for (size_t ll_c = 0; ll_c < m_NumOuterLayerCells; ll_c++) {
  //    //.... copy all content
  //    m_TransOuterCellsLayers[ll_c] = m_OuterCellsLayers[ll_c];
  //    //.... replace data pointer
  //    if(m_OuterCellsLayers[ll_c].m_Data == pivot_patch->getData()) {
  //      m_TransOuterCellsLayers[ll_c].m_Data = pivot_patch->getGpuData();
  //    } else {
  //      bool found = false;
  //      for(size_t i_p = 0; i_p < num_patches; i_p++) {
  //        Patch* patch = m_PatchGrid->getPatch(i_p);
  //        if(m_OuterCellsLayers[ll_c].m_Data == patch->getData()) {
  //          m_TransOuterCellsLayers[ll_c].m_Data = patch->getData();
  //          pivot_patch = patch;
  //          found = true;
  //          break;
  //        }
  //      }
  //      if(!found) {BUG;}
  //    }
  //  }

  cudaMalloc(&m_GpuInnerCellsLayers, sizeof(LSLayerDataExtrapol)*m_NumInnerLayerCells);
  CUDA_CHECK_ERROR;
  cudaMemcpy(m_GpuInnerCellsLayers, m_TransInnerCellsLayers, sizeof(LSLayerDataExtrapol)*m_NumInnerLayerCells, cudaMemcpyHostToDevice);
  CUDA_CHECK_ERROR;
  delete[] m_TransInnerCellsLayers;

  //  cudaMalloc(&m_GpuOuterCellsLayers, sizeof(LSLayerDataExtrapol)*m_NumOuterLayerCells);
  //  CUDA_CHECK_ERROR;
  //  cudaMemcpy(m_GpuOuterCellsLayers, m_TransOuterCellsLayers, sizeof(LSLayerDataExtrapol)*m_NumOuterLayerCells, cudaMemcpyHostToDevice);
  //  CUDA_CHECK_ERROR;
  //  delete[] m_TransOuterCellsLayers;
}


template <typename OP>
void GPU_LevelSetObjectBC<OP>::operator()()
{
  real relax = 0.5;

  ???
  GPU_LevelSetObjectBC_access_kernel(m_Op,
                                     m_Field,
                                     m_AbuseField,
                                     ???m_GpuInnerCellsLayers???);

  // hier müssen alle mit GPU_LevelSetObjectBC_access_kernel durch sein

  ???
  GPU_LevelSetObjectBC_operateInner_kernel(m_Op,
                                           m_Field,
                                           m_AbuseField,
                                           relax,
                                           ???m_GpuInnerCellsLayers???);

}

#endif // GPU_LEVELSETOBJECTBC_H
