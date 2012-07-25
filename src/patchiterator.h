#ifndef PATCHITERATOR_H
#define PATCHITERATOR_H

#include "patch.h"

class PatchIterator
{

private: // attributes

  Patch* m_Patch;

public:

  PatchIterator(Patch *patch);

  Patch* patch() { return m_Patch; }

  virtual void compute(real factor) = 0;

};

#endif // PATCHITERATOR_H
