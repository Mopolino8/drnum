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
#ifndef UPWIND1_H
#define UPWIND1_H

#include "cartesianpatch.h"

struct Upwind1
{
  template <typename PATCH>
  CUDA_DH void project(PATCH *patch, real* var, size_t i_field, size_t num_vars,
                       size_t i1, size_t j1, size_t k1,
                       size_t,    size_t,    size_t,
                       real, real, real,
                       real, real, real)
  {
    for (size_t i_var = 0; i_var < num_vars; ++i_var) {
      var[i_var] = patch->f(i_field, i_var, i1, j1, k1);
    }
  }
};

#endif // UPWIND1_H
