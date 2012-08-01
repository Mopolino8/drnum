#include "patchiterator.h"

PatchIterator::PatchIterator(Patch *patch)
{
  m_Patch = patch;
}

void PatchIterator::transformShapes()
{
  for (list<Shape*>::iterator i = m_Shapes.begin(); i != m_Shapes.end(); ++i) {
    (*i)->reset();
    (*i)->transform(m_Patch->getTransformation());
  }
}
