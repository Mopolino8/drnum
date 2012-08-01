#ifndef SHAPE_H
#define SHAPE_H

#include "blockcfd.h"
#include "transformation.h"

class Shape
{

public: // methods

  virtual void transform(const Transformation& transformation) = 0;
  virtual void reset() = 0;

};

#endif // SHAPE_H
