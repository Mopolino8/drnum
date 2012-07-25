#ifndef TIMEINTEGRATION_H
#define TIMEINTEGRATION_H

#include "patchiterator.h"
#include <list>

class TimeIntegration
{

private: // attributes

  list<PatchIterator*> m_Iterators;


protected:

  void copyField(size_t i_src, size_t i_dst);
  void computeIterators(real factor);


public:

  TimeIntegration();

  void addIterator(PatchIterator *patch_iterator) { m_Iterators.push_back(patch_iterator); }

  virtual void operator()(real dt) = 0;

};

#endif // TIMEINTEGRATION_H
