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
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef MUTEX_F_H
#define MUTEX_F_H

#include "mutex.h"

FORTRAN(mutex_new) (int *i, int *id, int *is_onwer)
{
  *i = new_pointer();
  try {
    global_pointer[*i] = new Mutex(*id, *is_onwer);
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(mutex_delete) (int *i)
{
  try {
    delete global_pointer[*i];
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(mutex_lock) (int *i)
{
  try {
    Mutex *mutex;
    get_pointer(*i, mutex);
    mutex->lock();
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(mutex_unlock) (int *i)
{
  try {
    Mutex *mutex;
    get_pointer(*i, mutex);
    mutex->unlock();
  } catch (IpcException e) {
    e.print();
  }
}


#endif // MUTEX_F_H
