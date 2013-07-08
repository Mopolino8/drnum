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
