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
#include <iostream>

#include "sharedmemory.h"

using namespace std;

int main()
{
  SharedMemory shm(1);
  Mutex mutex(1);
  double *real_data;
  int *int_data;
  mutex.lock();
  shm.readArray("real_data", real_data);
  shm.readArray("int_data", int_data);
  double real_sum = 0;
  int int_sum = 0;
  for (int i = 0; i < 20; ++i) {
    real_sum += real_data[i];
    int_sum += int_data[i];
  }
  cout << "real_sum = " << real_sum << endl;
  cout << "int_sum  = " << int_sum << endl;
  string text = shm.readString("text1");
  cout << "text     = '" << text << "'" << endl;
  shm.readValue("real_sum", real_sum);
  shm.readValue("int_sum", int_sum);
  cout << "real_sum = " << real_sum << endl;
  cout << "int_sum  = " << int_sum << endl;
  mutex.unlock();
  int int_value = 1;
  do {
    int i;
    mutex.lock();
    shm.readValue("int_value", i);
    mutex.unlock();
    if (i != int_value) {
      cout << "new int_value = " << i << endl;
      int_value = i;
    }
  } while (int_value);
}
