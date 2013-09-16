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
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef BARRIER_F_H
#define BARRIER_F_H

#include "barrier.h"

FORTRAN(barrier_new) (int *i, int *id, int *is_owner)
{
  *i = new_pointer();
  try {
    global_pointer[*i] = new Barrier(*id, *is_owner);
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(barrier_delete) (int *i)
{
  try {
    delete global_pointer[*i];
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(barrier_wait) (int *i)
{
  try {
    Barrier *barrier;
    get_pointer(*i, barrier);
    barrier->wait();
  } catch (IpcException e) {
    e.print();
  }
}

#endif // BARRIER_F_H
