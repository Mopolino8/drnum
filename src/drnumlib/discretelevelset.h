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

#ifndef DISCRETELEVELSET_H
#define DISCRETELEVELSET_H

#include "drnum.h"

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>
//#include <CGAL/exceptions.h>

#include <vtkPolyData.h>
#include <vtkSTLReader.h>

#include "patchgrid.h"
#include "postprocessingvariables.h"

template <unsigned int DIM, unsigned int IVAR>
class DiscreteLevelSet
{

private: // types

  typedef CGAL::Simple_cartesian<double> K;

  typedef K::FT         FT;
  typedef K::Ray_3      Ray;
  typedef K::Line_3     Line;
  typedef K::Point_3    Point;
  typedef K::Segment_3  Segment;
  typedef K::Triangle_3 Triangle;

  typedef QVector<Triangle>::iterator                       TriangleIterator;
  typedef CGAL::AABB_triangle_primitive<K,TriangleIterator> TrianglePrimitive;
  typedef CGAL::AABB_traits<K, TrianglePrimitive>           TriangleTraits;
  typedef CGAL::AABB_tree<TriangleTraits>                   TriangleTree;
  typedef TriangleTree::Point_and_primitive_id              TrianglePointAndPrimitiveId;

  typedef QVector<Segment>::iterator                       SegmentIterator;
  typedef CGAL::AABB_segment_primitive<K, SegmentIterator> SegmentPrimitive;
  typedef CGAL::AABB_traits<K, SegmentPrimitive>           SegmentTraits;
  typedef CGAL::AABB_tree<SegmentTraits>                   SegmentTree;

  typedef TriangleTree::Intersection_and_primitive_id<Ray>::Type Intersection;


private: // attributes

  PatchGrid* m_PatchGrid;

  //SegmentTree          m_SegmentTree;
  //vtkPolyData*         m_PolyData;
  //QVector<vtkIdType>   m_Tri2Grid;
  //QVector<double>      m_Radius;      ///< Surface radius for mesh resolution.


protected: //

  void computeLevelSet(vtkPolyData* poly);


public:

  DiscreteLevelSet(PatchGrid *patch_grid);

  void readStlGeometry(QString stl_file_name);

};


template <unsigned int IVAR>
class LevelSetPlotVars : public PostProcessingVariables
{
public:

  virtual int numScalars() const { return 1; }
  virtual int numVectors() const { return 0; }

  virtual string getScalarName(int) const { return "G"; }
  virtual string getVectorName(int) const { BUG; }
  virtual real   getScalar(int, real* var) const { return var[IVAR]; }
  virtual vec3_t getVector(int, real*)     const { BUG; return vec3_t(0,0,0); }
};


template <unsigned int DIM, unsigned int IVAR>
DiscreteLevelSet<DIM,IVAR>::DiscreteLevelSet(PatchGrid *patch_grid)
{
  m_PatchGrid = patch_grid;
}

template <unsigned int DIM, unsigned int IVAR>
void DiscreteLevelSet<DIM,IVAR>::readStlGeometry(QString stl_file_name)
{
  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName(qPrintable(stl_file_name));
  reader->MergingOff();
  reader->Update();
  vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
  poly->DeepCopy(reader->GetOutput());
  poly->BuildCells();
  computeLevelSet(poly);
}

template <unsigned int DIM, unsigned int IVAR>
void DiscreteLevelSet<DIM,IVAR>::computeLevelSet(vtkPolyData *poly)
{
  QVector<Triangle> triangles;
  QVector<vec3_t>   normals;
  TriangleTree      triangle_tree;
  //QVector<Segment>     m_Segments;

  // build triangle tree
  {
    int num_faces = poly->GetNumberOfCells();
    triangles.clear();
    normals.clear();
    triangles.fill(Triangle(), num_faces);
    normals.fill(vec3_t(0,0,0), num_faces);
    for (vtkIdType id_cell = 0; id_cell < num_faces; ++id_cell) {
      vtkIdType num_pts, *pts;
      poly->GetCellPoints(id_cell, num_pts, pts);
      if (num_pts != 3) {
        ERROR("only triangulated geometries are allowed");
      }
      dvec3_t a, b, c;
      poly->GetPoint(pts[0], a.data());
      poly->GetPoint(pts[1], b.data());
      poly->GetPoint(pts[2], c.data());
      triangles[id_cell] = Triangle(Point(a[0], a[1], a[2]), Point(b[0], b[1], b[2]), Point(c[0], c[1], c[2]));
      vec3_t u = b - a;
      vec3_t v = c - a;
      normals[id_cell] = u.cross(v);
      normals[id_cell].normalise();
    }
    triangle_tree.rebuild(triangles.begin(), triangles.end());
    triangle_tree.accelerate_distance_queries();
  }

  int count = 0;
  real var[DIM];
  dim_t<DIM> dim;
  for (size_t i_patch = 0; i_patch < m_PatchGrid->getNumPatches(); ++i_patch) {
    for (size_t i_cell = 0; i_cell < m_PatchGrid->getPatch(i_patch)->variableSize(); ++i_cell) {
      vec3_t x = m_PatchGrid->getPatch(i_patch)->xyzoCell(i_cell);
      try {
        Point p(x[0], x[1], x[2]);
        TrianglePointAndPrimitiveId result = triangle_tree.closest_point_and_primitive(p);
        Point cp = result.first;
        Triangle* T = result.second;
        int id = (T - triangles.begin());
        vec3_t x_snap = vec3_t(cp[0], cp[1], cp[2]);
        m_PatchGrid->getPatch(i_patch)->getVar(dim, 0, i_cell, var);
        vec3_t v = x - x_snap;
        real   g = v.abs();
        if (v*normals[id] < 0) {
          g *= -1;
        }
        var[IVAR] = g;
        m_PatchGrid->getPatch(i_patch)->setVar(dim, 0, i_cell, var);
      } catch (...) {
        ERROR("cannot comput distance");
      }
      ++count;
      if (count >= 10000) {
        cout << i_cell + 1 << endl;
        count = 0;
      }
    }
  }
}


#endif // DISCRETELEVELSET_H
