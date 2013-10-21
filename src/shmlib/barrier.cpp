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
#include "barrier.h"

Barrier::Barrier(int id_num, bool is_owner) : Semaphore(id_num, is_owner, 2)
{
  if (isOwner()) {
    set(0,1);
    set(1,1);
  } else {
    incr(0);
    incr(1);
  }
}

Barrier::~Barrier()
{
  wait();
  decr(0);
  decr(1);
}

void Barrier::wait()
{
  using namespace std;
  int s0 = get(0);
  int s1 = get(1);
  int max_count = get(1);
  decr(0);
  Semaphore::wait(0);
  decr(1);
  if (isOwner()) {
    Semaphore::wait(1);
    set(0, max_count);
    set(1, max_count);
  }
}

void Barrier::release()
{
  using namespace std;
  if (!isOwner() && get(0) != get(1)) {
    wait();
  }
}

void Barrier::print()
{
  std::cout << get(0) << " " << get(1) << std::endl;
}
