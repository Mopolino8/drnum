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
#ifndef GPU_RUNGEKUTTA_H
#define GPU_RUNGEKUTTA_H

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "rungekutta.h"
#include "iterators/gpu_patchiterator.h"

class GPU_RungeKutta : public RungeKutta
{

protected:

  virtual void copyField(size_t i_src, size_t i_dst)
  {
    for (list<PatchIterator*>::iterator i = m_Iterators.begin(); i != m_Iterators.end(); ++i) {
      GPU_PatchIterator* gpu_patch_iterator = dynamic_cast<GPU_PatchIterator*>(*i);
      if (gpu_patch_iterator == NULL) {
        BUG;
      }
      gpu_patch_iterator->copyField(i_src, i_dst);
    }
  }

};

#endif // GPU_RUNGEKUTTA_H
