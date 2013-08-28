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
#include "blockcfd.h"

unsigned long int global_flops     = 0;
unsigned long int global_flops_x86 = 0;

time_t global_start_time = time(NULL);
time_t global_last_time  = time(NULL);

void startTiming()
{
  global_flops     = 0;
  global_flops_x86 = 0;
  global_start_time = time(NULL);
  global_last_time = global_start_time;
  cout << "\nTiming started.\n" << endl;
}

void printTiming()
{
  time_t now = time(NULL);
  cout << now - global_start_time << " seconds since timing began" << endl;
  cout << now - global_last_time << " seconds since last timing output" << endl;
  global_last_time = now;
}

void stopTiming()
{
  int secs = time(NULL) - global_start_time;
  cout << "\nTiming finished:\n";
  if (global_flops > 0) {
    cout << "  floating point operations                : " << global_flops/1000000     << "*10^6\n";
    cout << "  floating point operations (X86 weighted) : " << global_flops_x86/1000000 << "*10^6\n";
    cout << "  floating point operations                : " << global_flops/1000000/secs     << " MFlops\n";
    cout << "  floating point operations (X86 weighted) : " << global_flops_x86/1000000/secs << " MFlops\n";
    cout << "  seconds                                  : " << secs << "\n";
  } else {
    cout << "  seconds : " << secs << "\n";
  }
  cout << endl;
}


