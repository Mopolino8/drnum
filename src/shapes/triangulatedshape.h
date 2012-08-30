#ifndef TRIANGULATEDSHAPE_H
#define TRIANGULATEDSHAPE_H

#include "shapes/shape.h"

#include <string>
#include <vtkCellLocator.h>

class TriangulatedShape : public Shape
{

private:

  vtkPolyData    *m_PolyDataCopy;
  vtkPolyData    *m_PolyData;
  vtkCellLocator *m_Locator;

public:

  TriangulatedShape();
  ~TriangulatedShape();

  void readStlFile(string file_name);

  bool getBoundaryMetric(real x1, real y1, real z1,
                         real x2, real y2, real z2,
                         real &k, real &nx, real &ny, real &nz);
  bool isInside(real x, real y, real z);

  virtual void transform(const Transformation &transformation);
  virtual void reset();

};

#endif // TRIANGULATEDSHAPE_H
