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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef TIMEINTEGRATION_H
#define TIMEINTEGRATION_H

#include "iterators/patchiterator.h"
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
  void copyDonorData(size_t i_field);

  virtual void operator()(real dt) = 0;

};

#endif // TIMEINTEGRATION_H
