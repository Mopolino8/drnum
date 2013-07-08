#include "timeintegration.h"

TimeIntegration::TimeIntegration()
{
}


void TimeIntegration::copyField(size_t i_src, size_t i_dst)
{
  for (list<PatchIterator*>::iterator i = m_Iterators.begin(); i != m_Iterators.end(); ++i) {
    (*i)->copyField(i_src, i_dst);
  }
}


void TimeIntegration::computeIterators(real factor)
{
  for (list<PatchIterator*>::iterator i = m_Iterators.begin(); i != m_Iterators.end(); ++i) {
    (*i)->computeAll(factor);
  }
}

void TimeIntegration::copyDonorData(size_t i_field)
{
  for (list<PatchIterator*>::iterator i = m_Iterators.begin(); i != m_Iterators.end(); ++i) {
    (*i)->copyDonorData(i_field);
  }
}

void TimeIntegration::runPostOperations()
{
  for (list<GenericOperation*>::iterator i = m_PostOperations.begin(); i != m_PostOperations.end(); ++i) {
    (*(*i))();
  }
}

