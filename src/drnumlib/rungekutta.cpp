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
#include "rungekutta.h"

RungeKutta::RungeKutta()
{
}

void RungeKutta::operator()(real dt)
{
  bool first_step = true;
  copyDonorData(0);
  copyField(0, 1);
  for (list<real>::iterator i = m_Alpha.begin(); i != m_Alpha.end(); ++i) {
    if (!first_step) {
      copyDonorData(0);
    }
    first_step = false;
    computeIterators((*i)*dt);
    runPostOperations();
    countFlops(1);
  }
}
