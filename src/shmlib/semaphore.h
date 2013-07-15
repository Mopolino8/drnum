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
#ifndef SEMAPHORE_H
#define SEMAPHORE_H

class Semaphore;

#include "ipcobject.h"

class Semaphore : public IPCObject
{

#ifndef NO_IPC
  union semun
  {
    int val;
    struct semid_ds *buf;
    unsigned short int *array;
    struct seminfo *__buf;
  };
#endif

  int  m_NumSems;

protected:

  void decr(int sem_num);
  void incr(int sem_num);
  void wait(int sem_num);
  int  get(int sem_num);
  void set(int sem_num, int value);

public:

  Semaphore(int id_num, bool is_owner = false, int num_sems = 1);
  ~Semaphore();

  virtual void close();

  int  numSems() { return m_NumSems; }

};

#endif // SEMAPHORE_H
