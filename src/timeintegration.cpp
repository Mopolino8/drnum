#include "timeintegration.h"

TimeIntegration::TimeIntegration()
{
}


void TimeIntegration::copyField(size_t i_src, size_t i_dst)
{
  for (list<PatchIterator*>::iterator i = m_Iterators.begin(); i != m_Iterators.end(); ++i) {
    (*i)->patch()->copyField(i_src, i_dst);
  }
}


void TimeIntegration::computeIterators(real factor)
{
  for (list<PatchIterator*>::iterator i = m_Iterators.begin(); i != m_Iterators.end(); ++i) {
    (*i)->compute(factor);
  }
}

void TimeIntegration::runPostOperations()
{
  for (list<GenericOperation*>::iterator i = m_PostOperations.begin(); i != m_PostOperations.end(); ++i) {
    (*(*i))();
  }
}

