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
// + DrNUM is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef COMPRESSIBLEVARIABLESANDG_H
#define COMPRESSIBLEVARIABLESANDG_H

#include "blockcfd.h"
#include "postprocessingvariables.h"

#include <QString>


/** @todo
  * This is an intermediate test variable set to allow post processing on blockobject
  * design and simulations.
  */

template <typename TGas>
class CompressibleVariablesAndG : public PostProcessingVariables
{

  TGas m_Gas;
  real m_RefPressure;
  real m_RefTemperature;

public:
  CompressibleVariablesAndG();
  virtual int numScalars() const { return 6; }
  virtual int numVectors() const { return 1; }

  void setReferenceTemperature(real T) { m_RefTemperature = T; }
  void setReferencePressure(real p)    { m_RefPressure = p; }

  virtual string getScalarName(int i) const;
  virtual string getVectorName(int i) const;
  virtual real   getScalar(int i, real* var) const;
  virtual vec3_t getVector(int i, real* var) const;

};

template <typename TGas>
CompressibleVariablesAndG<TGas>::CompressibleVariablesAndG()
{
  m_RefPressure = 1e5;
  m_RefTemperature = 300;
}

template <typename TGas>
string CompressibleVariablesAndG<TGas>::getScalarName(int i) const
{
  if (i == 0) return "Ma";
  if (i == 1) return "p";
  if (i == 2) return "T";
  if (i == 3) return "rho";
  if (i == 4) return "S";
  if (i == 5) return "G";   // grey value. Body contours: G = 0.5
  BUG;
  return "N/A";
}

template <typename TGas>
string CompressibleVariablesAndG<TGas>::getVectorName(int i) const
{
  if (i == 0) return "U";
  BUG;
  return "N/A";
}

template <typename TGas>
real CompressibleVariablesAndG<TGas>::getScalar(int i, real *var) const
{
  real p, T, u, v, w;
  m_Gas.conservativeToPrimitive(var, p, T, u, v, w);

  if (i == 0) return sqrt((u*u + v*v + w*w)/(m_Gas.gamma()*m_Gas.R()*T));
  if (i == 1) return p;
  if (i == 2) return T;
  if (i == 3) return var[0];
  if (i == 4) return TGas::cp(var)*log(T/m_RefTemperature) - TGas::R(var)*log(p/m_RefPressure);
  if (i == 5) return var[5];
  BUG;
  return 0;
}

template <typename TGas>
vec3_t CompressibleVariablesAndG<TGas>::getVector(int i, real *var) const
{
  real p, T, u, v, w;
  m_Gas.conservativeToPrimitive(var, p, T, u, v, w);

  if (i == 0) return vec3_t(u, v, w);
  BUG;
  return vec3_t(0,0,0);
}

#endif // COMPRESSIBLEVARIABLESANDG_H
