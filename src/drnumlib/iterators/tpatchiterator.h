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
