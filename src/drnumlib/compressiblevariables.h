#ifndef COMPRESSIBLEVARIABLES_H
#define COMPRESSIBLEVARIABLES_H

#include "blockcfd.h"
#include "postprocessingvariables.h"

#include <QString>

template <typename TGas>
class CompressibleVariables : public PostProcessingVariables
{

  TGas m_Gas;
  real m_RefPressure;
  real m_RefTemperature;

public:

  CompressibleVariables();

  virtual int numScalars() const { return 5; }
  virtual int numVectors() const { return 1; }

  void setReferenceTemperature(real T) { m_RefTemperature = T; }
  void setReferencePressure(real p)    { m_RefPressure = p; }

  virtual string getScalarName(int i) const;
  virtual string getVectorName(int i) const;
  virtual real   getScalar(int i, real* var) const;
  virtual vec3_t getVector(int i, real* var) const;

};

template <typename TGas>
CompressibleVariables<TGas>::CompressibleVariables()
{
  m_RefPressure = 1e5;
  m_RefTemperature = 300;
}

template <typename TGas>
string CompressibleVariables<TGas>::getScalarName(int i) const
{
  if (i == 0) return "Ma";
  if (i == 1) return "p";
  if (i == 2) return "T";
  if (i == 3) return "rho";
  if (i == 4) return "S";
  BUG;
  return "N/A";
}

template <typename TGas>
string CompressibleVariables<TGas>::getVectorName(int i) const
{
  if (i == 0) return "U";
  BUG;
  return "N/A";
}

template <typename TGas>
real CompressibleVariables<TGas>::getScalar(int i, real *var) const
{
  real p, T, u, v, w;
  m_Gas.conservativeToPrimitive(var, p, T, u, v, w);

  if (i == 0) return sqrt((u*u + v*v + w*w)/(m_Gas.gamma()*m_Gas.R()*T));
  if (i == 1) return p;
  if (i == 2) return T;
  if (i == 3) return var[0];
  if (i == 4) return TGas::cp(var)*log(T/m_RefTemperature) - TGas::R(var)*log(p/m_RefPressure);
  BUG;
  return 0;
}

template <typename TGas>
vec3_t CompressibleVariables<TGas>::getVector(int i, real *var) const
{
  real p, T, u, v, w;
  m_Gas.conservativeToPrimitive(var, p, T, u, v, w);

  if (i == 0) return vec3_t(u, v, w);
  BUG;
  return vec3_t(0,0,0);
}

#endif // COMPRESSIBLEVARIABLES_H
