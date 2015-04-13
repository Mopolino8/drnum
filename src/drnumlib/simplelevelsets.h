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

#ifndef SIMPLELEVELSETS_H
#define SIMPLELEVELSETS_H

#include "cartesianpatch.h"
#include "gpu_cartesianpatch.h"
#include "postprocessingvariables.h"

template <typename LS>
struct GenericLevelSetPlotVars : public PostProcessingVariables
{
public:

  virtual int numScalars() const { return 1; }
  virtual int numVectors() const { return 0; }

  virtual string getScalarName(int) const { return "G"; }
  virtual string getVectorName(int) const { BUG; }
  virtual real   getScalar(int, real*, vec3_t x) const { return LS::G(x[0], x[1], x[2]); }
  virtual vec3_t getVector(int, real*, vec3_t)   const { BUG; return vec3_t(0,0,0); }

};


template <int T_X0, int T_Y0, int T_Z0, int T_R, int T_H>
struct LevelSetXCylinder
{
  CUDA_DH static real G(real x, real y, real z)
  {
    real R = 1e-3*T_R;
    real H = 1e-3*T_H;
    x -= 1e-3*T_X0;
    y -= 1e-3*T_Y0;
    z -= 1e-3*T_Z0;
    real h = x;
    real r = sqrt(y*y + z*z);
    if (h >= 0 && h <= H) {
      if (r <= R) return -min(R-r, min(h, H-h));
      else        return r-R;
    } else {
      if (r <= R) {
        if (h < 0) return -h;
        else       return  h-H;
      } else {
        real dh = min(-h, h-H);
        real dr = R - r;
        return sqrt(dh*dh + dr*dr);
      }
    }
  }

  template <typename TPatch>
  CUDA_DH static real G(TPatch& patch, size_t i, size_t j, size_t k, size_t = 0)
  {
    real x, y, z;
    patch.xyzIJK(i, j, k, x, y, z);
    return G(x, y, x);
  }

  template <typename TPatch>
  CUDA_DH static void updateG(TPatch&, size_t, size_t, size_t, real, size_t = 0) {}
};

#endif // SIMPLELEVELSETS_H
