#ifndef TPATCHITERATOR_H
#define TPATCHITERATOR_H

#include <vector>

#include "patchgrid.h"
#include "patchiterator.h"

template <typename T, unsigned int DIM, typename OP>
class TPatchIterator : public PatchIterator
{

protected:

  vector<T*> m_Patches;
  OP         m_Op;

public:

  TPatchIterator(OP op);

  //virtual void addPatch(T* patch);
  virtual void addPatch(Patch* patch);

};


template <typename T, unsigned int DIM, typename OP>
TPatchIterator<T, DIM, OP>::TPatchIterator(OP op) : PatchIterator()
{
  //m_Patches.reserve(max(size_t(100), patch_grid.getNumPatches()));
  m_Op = op;
}

template <typename T, unsigned int DIM, typename OP>
void TPatchIterator<T, DIM, OP>::addPatch(Patch* patch)
{
  /// @todo currently ugly. Find concept to avoid storing m_Patches twice. See PatchIterator.h
  PatchIterator::addPatch(patch);

  m_Patches.push_back(dynamic_cast<T*>(patch));
  //m_Patches.push_back(patch);
}


#endif // TPATCHITERATOR_H
