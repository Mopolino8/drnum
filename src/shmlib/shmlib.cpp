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
#include <iostream>
#include <cstdlib>
#include <string>
#include <cstring>

using namespace std;

#include "ipcobject.h"

#define FORTRAN(NAME) extern "C" void NAME ## _


#define NUMBER_OF_POINTERS 1000

IPCObject* global_pointer[NUMBER_OF_POINTERS];
int last_pointer = 0;

template <class T>
void get_pointer(int i, T* &pointer)
{
  pointer = dynamic_cast<T*>(global_pointer[i]);
  if (!pointer) {
    cerr << "pointer conversion error: index=" << i << endl;
    exit(EXIT_FAILURE);
  }
}

int new_pointer()
{
  if (last_pointer >= NUMBER_OF_POINTERS) {
    cerr << "FORTRAN shared memory interface: 'number of created objects exceeded'" << endl;
    exit(EXIT_FAILURE);
  }
  ++last_pointer;
  return last_pointer - 1;
}

#include "mutex.f.h"
#include "barrier.f.h"
#include "sharedmemory.f.h"



FORTRAN(cpp_test) (char *text)
{
  int i = 0;
  string s = "";
  while (text[i] != 32 && i < 20) {
    cout << text[i] << ", " << int(text[i]) << endl;
    s += text[i];
    ++i;
  }
  cout << s << endl;
}

