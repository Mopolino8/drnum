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
