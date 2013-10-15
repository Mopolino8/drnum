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
#include "sharedmemory.h"

SharedMemory::SharedMemory(int id_num, int size, bool is_owner) : IPCObject(id_num, is_owner)
{
#ifndef NO_IPC
  if (m_Owner) {
    m_Size = size;
    setId(shmget(key(), size, IPC_CREAT | IPC_EXCL | 0660));
    if (id() == -1) {
      setId(shmget(key(), 1, IPC_CREAT | 0660));
      if (id() == -1) {
        throw IpcException("unable to create sshared memory\nA\nshmget " + errorText(errno));
      }
      SharedMemory::close();
      setId(shmget(key(), size, IPC_CREAT | IPC_EXCL | 0660));
      if (id() == -1) {
        throw IpcException("unable to create sshared memory\nB\nshmget " + errorText(errno));
      }
    }
  } else {
    setId(shmget(key(), size, IPC_CREAT | 0660));
    if (id() == -1) {
      throw IpcException("unable to create sshared memory\nC\nshmget " + errorText(errno));
    }
    struct shmid_ds buf;
    shmctl(id(), IPC_STAT, &buf);
    m_Size = buf.shm_segsz;
  }
  m_Buffer = (char*) shmat(id(), 0, 0);
  if (!m_Buffer) {
    throw IpcException("shmat " + errorText(errno));
  }
  m_MaxNameLength = 32;
  m_MaxNumArrays = 100;
  m_ArrayDescrLength = m_MaxNameLength + 3*sizeof(int);
  m_Offset = 2*sizeof(int);
  if (m_Owner) {
    reset();
  }
#else
  size = size;
#endif
}

SharedMemory::~SharedMemory()
{
  SharedMemory::close();
}

void SharedMemory::close()
{
#ifndef NO_IPC
  if (isOwner()) {
    shmctl(id(), 0, IPC_RMID);
  }
#endif
}

void SharedMemory::reset()
{
  int *array_length = 0;
  int *data_index = 0;
  get_pointer(0, array_length);
  get_pointer(4, data_index);
  *array_length = 0;
  *data_index = m_Offset + m_MaxNumArrays*m_ArrayDescrLength;
}

size_t SharedMemory::numArrays()
{
  int *N = 0;
  get_pointer(0, N);
  return *N;
}

size_t SharedMemory::dataIndex()
{
  int *I = 0;
  get_pointer(4, I);
  return *I;
}

std::string SharedMemory::arrayName(int i)
{
  char *C = 0;
  get_pointer(indexOfName(i), C);
  size_t L = 0;
  std::string name = "";
  while (*C != char(0) && L < m_MaxNameLength) {
    name += *C;
    ++C;
    ++L;
  }
  return name;
}

int SharedMemory::arrayIndex(std::string name)
{  
  if (name == "dn2of-data-0") {
    std::cout << "SharedMemory::arrayIndex(\"" << name << "\")" << std::endl;
  }
  size_t i = 0;
  while (i < numArrays() && arrayName(i) != name) {
    ++i;
  }
  if (i >= numArrays()) {
    return -1;
  }
  return i;
}

std::string SharedMemory::readString(std::string name)
{
  int i = arrayIndex(name);
  int L = arrayLength(i);
  char *array = new char[L];
  readArray(name, array);
  std::string s = "";
  --L;
  while (L > 0 && array[L] == ' ') {
    --L;
  }
  for (int j = 0; j <= L; ++j) {
    if (array[j] == char(0)) {
      break;
    }
    s += array[j];
  }
  delete [] array;
  return s;
}
