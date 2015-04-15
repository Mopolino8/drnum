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
#include "postprocessingvariables.h"

#ifdef GPU
#include "gpu_cartesianpatch.h"
#endif

template <typename LS>
struct GenericLevelSetPlotVars : public PostProcessingVariables
{

protected:

  LS& m_Ls;


public:

  GenericLevelSetPlotVars(LS& ls) : PostProcessingVariables(), m_Ls(ls) {}

  virtual int numScalars() const { return 1; }
  virtual int numVectors() const { return 0; }

  virtual string getScalarName(int) const { return "G"; }
  virtual string getVectorName(int) const { BUG; }
  virtual real   getScalar(int, real*, vec3_t x) const { return m_Ls.G(x[0], x[1], x[2]); }
  virtual vec3_t getVector(int, real*, vec3_t)   const { BUG; return vec3_t(0,0,0); }

};


class LevelSetXCylinder
{

protected: // attributes

  real m_X0;
  real m_Y0;
  real m_Z0;
  real m_R;
  real m_H;


public: // methods

  CUDA_DH LevelSetXCylinder() {}

  CUDA_DH LevelSetXCylinder(real x0, real y0, real z0, real R, real H)
  {
    m_X0 = x0;
    m_Y0 = y0;
    m_Z0 = z0;
    m_R = R;
    m_H = H;
  }

  CUDA_DH real G(real x, real y, real z)
  {
    x -= m_X0;
    y -= m_Y0;
    z -= m_Z0;
    real h = x;
    real r = sqrt(y*y + z*z);
    if (h >= 0 && h <= m_H) {
      if (r <= m_R) return -min(m_R-r, min(h, m_H-h));
      else          return r-m_R;
    } else {
      if (r <= m_R) {
        if (h < 0) return -h;
        else       return  h-m_H;
      } else {
        real dh = min(-h, h-m_H);
        real dr = m_R - r;
        return sqrt(dh*dh + dr*dr);
      }
    }
  }

  template <typename TPatch>
  CUDA_DH real G(TPatch& patch, size_t i, size_t j, size_t k, size_t = 0)
  {
    real x, y, z;
    patch.xyzoIJK(i, j, k, x, y, z);
    return G(x, y, x);
  }

};

#endif // SIMPLELEVELSETS_H
