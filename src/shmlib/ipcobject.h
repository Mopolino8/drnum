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
#ifndef IPCOBJECT_H
#define IPCOBJECT_H

class IPCObject;

#include <sys/ipc.h>
#include <sys/sem.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <errno.h>
#include <string>
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>

#include "ipcexception.h"

class IPCObject
{

private:

  key_t m_Key;
  int   m_Id;

protected:

  bool    m_Owner;
  void    setId(int id) { m_Id = id; }

public:

  IPCObject(int id_num, bool is_owner = false);
  virtual ~IPCObject() {}

  void error(std::string msg);
  virtual void close() = 0;

  key_t key()     { return m_Key; }
  int   id()      { return m_Id; }
  bool  isOwner() { return m_Owner; }

  std::string errorText(int err);

};

#endif // IPCOBJECT_H
