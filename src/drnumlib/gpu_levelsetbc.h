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

#ifndef GPU_LEVELSETBC_H
#define GPU_LEVELSETBC_H

template <unsigned int DIM, typename TCPU, typename TGPU, typename BC>
class GPU_LevelSetBC;

#include "cudatools.h"

#include "drnum.h"
#include "genericoperation.h"

template <unsigned int DIM, typename TCPU, typename TGPU, typename BC>
class GPU_LevelSetBC : public GenericOperation
{

protected: // attributes

  PatchGrid*     m_PatchGrid;
  vector<TCPU*>  m_Patches;
  vector<TGPU>   m_GpuPatches;
  BC             m_Bc;


public: // methods

  GPU_LevelSetBC(PatchGrid* patch_grid);

};


template <unsigned int DIM, typename BC>
GPU_LevelSetBC<DIM,BC>::GPU_LevelSetBC(PatchGrid* patch_grid)
{
  m_PatchGrid = patch_grid;
  for (size_t i = 0; i < m_PatchGrid->getNumPatches(); ++i_patch) {
    TCPU* patch = dynamic_cast<TCPU*>(m_PatchGrid->getPatch(i));
    if (patch) {
      m_Patches.push_back(patch);
      m_GpuPatches.push_back(TGPU(patch));
    }
  }
}

#endif // GPU_LEVELSETBC_H
