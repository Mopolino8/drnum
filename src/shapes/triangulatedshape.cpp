#include "triangulatedshape.h"

#include <vtkSTLReader.h>

TriangulatedShape::TriangulatedShape()
{
  m_Locator = vtkCellLocator::New();
  m_PolyData = vtkPolyData::New();
  m_PolyDataCopy = vtkPolyData::New();
}

TriangulatedShape::~TriangulatedShape()
{
  m_Locator->Delete();
  m_PolyData->Delete();
  m_PolyDataCopy->Delete();
}

void TriangulatedShape::readStlFile(string file_name)
{
  vtkSTLReader* stl = vtkSTLReader::New();
  stl->SetFileName(file_name.c_str());
  stl->MergingOn();
  stl->Update();
  m_PolyData->DeepCopy(stl->GetOutput());
  m_PolyData->BuildCells();
  int nc = m_PolyData->GetNumberOfCells();
  int nn = m_PolyData->GetNumberOfPoints();
  m_PolyDataCopy->DeepCopy(stl->GetOutput());
  stl->Delete();
  m_Locator->SetDataSet(m_PolyData);
  m_Locator->BuildLocator();
}

bool TriangulatedShape::getBoundaryMetric(real x1, real y1, real z1,
                                          real x2, real y2, real z2,
                                          real &k,
                                          real &nx, real &ny, real &nz)
{
  double p1[3] = {x1, y1, z1};
  double p2[3] = {x2, y2, z2};
  double x[3], r[3];
  double t = 0;
  int id_sub = 0;
  vtkIdType id_cell = 0;
  if (!m_Locator->IntersectWithLine(p1, p2, 0, t, x, r, id_sub, id_cell)) {
    return false;
  }
  k = t;
  vtkIdType num_pts, *pts;
  m_PolyData->GetCellPoints(id_cell, num_pts, pts);
  if (num_pts != 3) {
    BUG;
  }
  dpvec3_t a, b, c;
  m_PolyData->GetPoint(pts[0], a.data());
  m_PolyData->GetPoint(pts[1], b.data());
  m_PolyData->GetPoint(pts[2], c.data());
  vec3_t u = b - a;
  vec3_t v = c - a;
  vec3_t n = u.cross(v);
  n.normalise();
  nx = n[0];
  ny = n[1];
  nz = n[2];
  return true;
}

bool TriangulatedShape::isInside(real x, real y, real z)
{
  static const vec3_t Dx(1,0,0);
  vec3_t n;
  real k;
  if (x > 1 && z < 0.2) {
    int a = 0;
  }
  if (!getBoundaryMetric(x, y, z, x + 100, y, z, k, n[0], n[1], n[2])) {
    return false;
  }
  if (n*Dx > 0) {
    return true;
  }
  return false;
}

void TriangulatedShape::transform(const Transformation &transformation)
{
  m_PolyDataCopy->DeepCopy(m_PolyData);
  for (vtkIdType id_node = 0; id_node < m_PolyData->GetNumberOfPoints(); ++id_node) {
    dpvec3_t x;
    m_PolyData->GetPoint(id_node, x.data());
    vec3_t xr(x[0], x[1], x[2]);
    xr = transformation(xr);
    x[0] = xr[0];
    x[1] = xr[1];
    x[2] = xr[2];
    m_PolyData->GetPoints()->SetPoint(id_node, x.data());
  }
  m_Locator->BuildLocator();
}

void TriangulatedShape::reset()
{
  m_Locator->SetDataSet(m_PolyData);
  m_PolyData->DeepCopy(m_PolyDataCopy);
}
