#ifndef PATCHITERATOR_H
#define PATCHITERATOR_H

#include "patch.h"
#include "patchgrid.h"
#include "shapes/shape.h"
#include "codestring.h"

class PatchIterator
{

private: // attributes

  vector<Patch*> m_Patches;

  CodeString m_SolverCodes;


public:

  PatchIterator();

  size_t numPatches() { return m_Patches.size(); }
  Patch* getPatch(size_t i) { return m_Patches[i]; }
  void   computeAll(real factor);

  virtual void addPatch(Patch* patch) { m_Patches.push_back(patch); }

  virtual void compute(real factor, const vector<size_t>& patches) = 0;
  virtual void copyField(size_t i_src, size_t i_dst);

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
  for (int i = 0; i < numPatches(); ++i) {
    patches[i] = i;
  }
  compute(factor, patches);
}

inline void PatchIterator::copyField(size_t i_src, size_t i_dst)
{
  for (int i = 0; i < numPatches(); ++i) {
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

#endif // PATCHITERATOR_H
