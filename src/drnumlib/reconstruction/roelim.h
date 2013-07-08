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
#ifndef ROELIM_H
#define ROELIM_H

#include "blockcfd.h"

struct RoeLim
{
  static CUDA_DH real delta(real delta1, real delta2)
  {
    countFlops(1);
    if (delta1*delta2 <= 0) {
      return 0;
    } else {
      countFlops(7);
      static const real eps = 1e-6;
      real D1  = nonZero(delta1, eps);
      real D2  = nonZero(delta2, eps);
      real lim = min(real(1), min(fabs(2*D1/(D1+D2)), fabs(2*D2/(D1+D2))));
      return lim*delta1;
    }
  }
};

#endif // ROELIM_H
