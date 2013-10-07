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
#include "ipcobject.h"

#include <sstream>

IPCObject::IPCObject(int id_num, bool is_owner)
{
  m_Key = 0;
#ifndef NO_IPC
  m_Key = ftok(".", id_num);
#endif
  m_Owner = is_owner;
}

#define trErr(E) if (err == E) return #E

std::string IPCObject::errorText(int err)
{
#ifndef NO_IPC
  trErr(EACCES);
  trErr(EEXIST);
  trErr(EIDRM);
  trErr(ENOENT);
  trErr(ENOMEM);
  trErr(ENOSPC);
  trErr(EINVAL);
  trErr(ERANGE);
  trErr(EPERM);
  trErr(EFAULT);

#endif
  std::stringstream out;
  out << err;
  return std::string("unknown error: ") + out.str();
}

#undef trErr

void IPCObject::error(std::string msg)
{
  close();
  throw IpcException(msg);
}

