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
#ifndef SHAREDMEMORY_F_H
#define SHAREDMEMORY_F_H

#include "sharedmemory.h"

std::string fortran2string(char *c)
{
  string s = "";
  while (*c != '.' && *c != 0x20) {
    s += *c;
    ++c;
  }
  return s;
}

FORTRAN(shm_new) (int *i, int *id, int *size, int *is_owner)
{
  *i = new_pointer();
  try {
    global_pointer[*i] = new SharedMemory(*id, *size, *is_owner);
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(shm_delete) (int *i)
{
  try {
    delete global_pointer[*i];
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(shm_reset) (int *i)
{
  try {
    SharedMemory *shm;
    get_pointer(*i, shm);
    shm->reset();
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(shm_write_real_array) (int *i, char *name, int *length, double *array)
{
  try {
    SharedMemory *shm;
    get_pointer(*i, shm);
    shm->writeArray(fortran2string(name), *length, array);
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(shm_write_integer_array) (int *i, char *name, int *length, int *array)
{
  try {
    SharedMemory *shm;
    get_pointer(*i, shm);
    shm->writeArray(fortran2string(name), *length, array);
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(shm_write_character_array) (int *i, char *name, int *length, char *array)
{
  try {
    SharedMemory *shm;
    get_pointer(*i, shm);
    shm->writeArray(fortran2string(name), *length, array);
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(shm_write_real_value) (int *i, char *name, double *value)
{
  try {
    SharedMemory *shm;
    get_pointer(*i, shm);
    shm->writeValue(fortran2string(name), value);
  } catch (IpcException e) {
    e.print();
  }
}

FORTRAN(shm_write_integer_value) (int *i, char *name, int *value)
{
  try {
    SharedMemory *shm;
    get_pointer(*i, shm);
    shm->writeValue(fortran2string(name), value);
  } catch (IpcException e) {
    e.print();
  }
}

#endif // SHAREDMEMORY_F_H
