#ifndef NOSHAPE_H
#define NOSHAPE_H

#include "shape.h"

class NoShape : public Shape
{

public: // methods

  bool getBoundaryMetric(real, real, real, real, real, real, real&, real&, real&, real&) { return false; }
  bool isInside(real x, real y, real z) { return false; }

  virtual void transform(const Transformation&) {}
  virtual void reset() {}

};

#endif // NOSHAPE_H
