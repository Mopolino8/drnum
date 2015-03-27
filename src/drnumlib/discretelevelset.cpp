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

#include "discretelevelset.h"

#include <CGAL/exceptions.h>
#include <vtkSTLReader.h>

DiscreteLevelSet::DiscreteLevelSet(PatchGrid *patch_grid, int var_index, QString stl_file_name)
{
  m_PatchGrid = patch_grid;
  m_VarIndex = var_index;
  vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
  reader->SetFileName(qPrintable(stl_file_name));
  reader->MergingOff();
  reader->Update();
  vtkSmartPointer<vtkPolyData> poly = vtkSmartPointer<vtkPolyData>::New();
  poly->DeepCopy(reader->GetOutput());
  poly->BuildCells();

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
  for (size_t i_patch = 0; i_patch < patch_grid->getNumPatches(); ++i_patch) {
    for (size_t i_cell = 0; i_cell < patch_grid->getPatch(i_patch)->variableSize(); ++i_cell) {
      vec3_t x = patch_grid->getPatch(i_patch)->xyzoCell(i_cell);
      try {
        Point p(x[0], x[1], x[2]);
        TrianglePointAndPrimitiveId result = triangle_tree.closest_point_and_primitive(p);
        Point cp = result.first;
        Triangle* T = result.second;
        int id = (T - triangles.begin());
        vec3_t x_snap = vec3_t(cp[0], cp[1], cp[2]);
      } catch (CGAL::Failure_exception) {
        ERROR("cannot comput distance");
      }
      ++count;
      if (count >= 10000) {
        cout << i_cell << endl;
        count = 0;
      }
    }
  }
}
