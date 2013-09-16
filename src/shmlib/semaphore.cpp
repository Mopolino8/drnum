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
#include "semaphore.h"

Semaphore::Semaphore(int id_num, bool is_owner, int num_sems) : IPCObject(id_num, is_owner)
{
  m_NumSems = num_sems;
#ifndef NO_IPC
  if (m_Owner) {
    setId(semget(key(), num_sems, IPC_CREAT | IPC_EXCL | 0660));
    if (id() == -1) {
      setId(semget(key(), num_sems, IPC_CREAT | 0660));
      if (id() == -1) {
        std::cout << "+++ 1 +++" << std::endl;
        std::cout << key() << ", " << num_sems << std::endl;
        throw IpcException("unable to create semaphore\nsemget " + errorText(errno));
      }
      Semaphore::close();
      setId(semget(key(), num_sems, IPC_CREAT | IPC_EXCL | 0660));
      if (id() == -1) {
        std::cout << "+++ 2 +++" << std::endl;
        throw IpcException("unable to create semaphore\nsemget " + errorText(errno));
        exit(EXIT_FAILURE);
      }
    }
    for (int i = 0; i < num_sems; ++i) {
      set(i, 1);
    }
  } else {
    setId(semget(key(), num_sems, IPC_CREAT | 0660));
    if (id() == -1) {
      std::cout << "+++ 3 +++" << std::endl;
      throw IpcException("unable to create semaphore\nsemget " + errorText(errno));
    }
  }
#endif
}

Semaphore::~Semaphore()
{
  Semaphore::close();
}

void Semaphore::close()
{
#ifndef NO_IPC
  if (isOwner()) {
    semctl(id(), 0, IPC_RMID, 0);
  }
#endif
}

void Semaphore::decr(int sem_num)
{
#ifndef NO_IPC
  struct sembuf sem_decr = {sem_num, -1, 0};
  semop(id(), &sem_decr, 1);
#endif
}

void Semaphore::incr(int sem_num)
{
#ifndef NO_IPC
  struct sembuf sem_incr_nw = {sem_num, 1, IPC_NOWAIT};
  semop(id(), &sem_incr_nw, 1);
#endif
}

void Semaphore::wait(int sem_num)
{
#ifndef NO_IPC
  struct sembuf sem_decr = {sem_num, 0, 0};
  semop(id(), &sem_decr, 1);
#endif
}

int Semaphore::get(int sem_num)
{
#ifndef NO_IPC
  return semctl(id(), sem_num, GETVAL, 0);
#endif
}

void Semaphore::set(int sem_num, int value)
{
#ifndef NO_IPC
  union semun semopts;
  semopts.val = value;
  if (semctl(id(), sem_num, SETVAL, semopts) == -1) {
    error("semctl " + errorText(errno));
  }
#endif
}

