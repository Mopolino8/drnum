#ifndef POSTPROCESSINGVARIABLES_H
#define POSTPROCESSINGVARIABLES_H

#include "blockcfd.h"

class PostProcessingVariables
{

public: // methods

  virtual int numScalars() const = 0;
  virtual int numVectors() const = 0;

  virtual string getScalarName(int i) const = 0;
  virtual string getVectorName(int i) const = 0;
  virtual real   getScalar(int i, real* var) const = 0;
  virtual vec3_t getVector(int i, real* var) const = 0;


};

#endif // POSTPROCESSINGVARIABLES_H
