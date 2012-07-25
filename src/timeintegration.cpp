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

