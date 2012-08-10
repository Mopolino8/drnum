#ifndef COMPRESSIBLEVARIABLES_H
#define COMPRESSIBLEVARIABLES_H

#include "blockcfd.h"

#include <QString>

template <typename TGas>
class CompressibleVariables
{

  TGas m_Gas;

public:

  int numScalars() const { return 3; }
  int numVectors() const { return 1; }

  QString getScalarName(int i) const;
  QString getVectorName(int i) const;
  real   getScalar(int i, real* var) const;
  vec3_t getVector(int i, real* var) const;

};

template <typename TGas>
QString CompressibleVariables<TGas>::getScalarName(int i) const
{
  if (i == 0) return "Ma";
  if (i == 1) return "p";
  if (i == 2) return "T";
  BUG;
}

template <typename TGas>
QString CompressibleVariables<TGas>::getVectorName(int i) const
{
  if (i == 0) return "U";
  BUG;
}

template <typename TGas>
real CompressibleVariables<TGas>::getScalar(int i, real *var) const
{
  real p, T, u, v, w;
  m_Gas.conservativeToPrimitive(var, p, T, u, v, w);

  if (i == 0) return sqrt((u*u + v*v + w*w)/(m_Gas.gamma()*m_Gas.R()*T));
  if (i == 1) return p;
  if (i == 2) return T;
  BUG;
}

template <typename TGas>
vec3_t CompressibleVariables<TGas>::getVector(int i, real *var) const
{
  real p, T, u, v, w;
  m_Gas.conservativeToPrimitive(var, p, T, u, v, w);

  if (i == 0) return vec3_t(u, v, w);
  BUG;
}

#endif // COMPRESSIBLEVARIABLES_H
