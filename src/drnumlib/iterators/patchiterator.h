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
#ifndef PATCHITERATOR_H
#define PATCHITERATOR_H

#include "patch.h"
#include "patchgrid.h"
#include "codestring.h"

class PatchIterator
{

private: // attributes

  vector<Patch*> m_Patches;
  vector<bool>   m_PatchActive;

  CodeString m_SolverCodes;


public:

  PatchIterator();

  virtual void updateHost() { BUG; }
  virtual void updateDevice() { BUG; }

  size_t numPatches() { return m_Patches.size(); }
  Patch* getPatch(size_t i) { return m_Patches[i]; }
  void   computeAll(real factor);

  virtual void addPatch(Patch* patch);

  virtual void compute(real factor, const vector<size_t>& patches) = 0;
  virtual void copyField(size_t i_src, size_t i_dst);
  virtual void copyDonorData(size_t i_field);

  void activatePatch(size_t i_patch);
  void deactivatePatch(size_t i_patch);
  bool patchActive(size_t i_patch) { return m_PatchActive[i_patch]; }

  /**
    * Set a CodeString as operation identifier.
    * @param code_string the CodeString to set.
    */
  void setCodeString(const CodeString& code_string);

  /**
    * Get
    * @return CodeString identifying operations.
    */
  CodeString getCodeString();

};

inline PatchIterator::PatchIterator()
{
  //m_Patches.reserve(max(size_t(100), patch_grid.getNumPatches()));
  m_SolverCodes = string("void");
}

inline void PatchIterator::computeAll(real factor)
{
  vector<size_t> patches(numPatches());
  for (size_t i = 0; i < numPatches(); ++i) {
    patches[i] = i;
  }
  compute(factor, patches);
}

inline void PatchIterator::copyField(size_t i_src, size_t i_dst)
{
  for (size_t i = 0; i < numPatches(); ++i) {
    m_Patches[i]->copyField(i_src, i_dst);
  }
}

inline void PatchIterator::setCodeString(const CodeString& code_string)
{
  m_SolverCodes = code_string;
}

inline CodeString PatchIterator::getCodeString()
{
  return m_SolverCodes;
}

inline void PatchIterator::copyDonorData(size_t i_field)
{
  for (size_t i_patch = 0; i_patch < m_Patches.size(); ++i_patch) {
    m_Patches[i_patch]->accessDonorDataDirect(i_field);
  }
}

inline void PatchIterator::addPatch(Patch *patch)
{
  m_Patches.push_back(patch);
  m_PatchActive.push_back(true);
}

inline void PatchIterator::activatePatch(size_t i_patch)
{
  if (i_patch >= m_Patches.size()) {
    BUG;
  }
  m_PatchActive[i_patch] = true;
}

inline void PatchIterator::deactivatePatch(size_t i_patch)
{
  if (i_patch >= m_Patches.size()) {
    BUG;
  }
  m_PatchActive[i_patch] = false;
}

#endif // PATCHITERATOR_H
