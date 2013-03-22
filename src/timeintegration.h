#ifndef TIMEINTEGRATION_H
#define TIMEINTEGRATION_H

#include "patchiterator.h"
#include "genericoperation.h"
#include <list>

class TimeIntegration
{

protected:

  list<PatchIterator*>    m_Iterators;
  list<GenericOperation*> m_PostOperations;

  virtual void copyField(size_t i_src, size_t i_dst);

  void computeIterators(real factor);
  void runPostOperations();


public:

  TimeIntegration();

  void addIterator(PatchIterator *patch_iterator) { m_Iterators.push_back(patch_iterator); }
  void addPostOperation(GenericOperation *operation) { m_PostOperations.push_back(operation); }

  virtual void operator()(real dt) = 0;

};

#endif // TIMEINTEGRATION_H
