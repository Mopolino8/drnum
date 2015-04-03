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

#include <vtkPolyData.h>
#include <vtkSTLReader.h>

#include "patchgrid.h"
#include "postprocessingvariables.h"
#include "cartesianpatch.h"
#include "gpu_cartesianpatch.h"

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

  PatchGrid*        m_PatchGrid;
  QVector<Triangle> m_Triangles;
  QVector<vec3_t>   m_Normals;
  TriangleTree      m_TriangleTree;
  real              m_Tol;

  //SegmentTree          m_SegmentTree;
  //vtkPolyData*         m_PolyData;
  //QVector<vtkIdType>   m_Tri2Grid;
  //QVector<double>      m_Radius;      ///< Surface radius for mesh resolution.


protected: //

  void computeLevelSet(vtkPolyData* poly);
  real computePointLevelSet(vec3_t x);
  void levelSetPerCell(size_t i_patch, int& count, int& total_count);
  void recursiveLevelSet(CartesianPatch* patch, size_t i_o, size_t j_o, size_t k_o, size_t i_m, size_t j_m, size_t k_m);
  void interpolate(CartesianPatch* patch, size_t i_o, size_t j_o, size_t k_o, size_t i_m, size_t j_m, size_t k_m);
  real distance(vec3_t x1, vec3_t x2);


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
  m_Tol = 0.55;
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
  //QVector<Segment>     m_Segments;

  // build triangle tree
  {
    int num_faces = poly->GetNumberOfCells();
    m_Triangles.clear();
    m_Normals.clear();
    m_Triangles.fill(Triangle(), num_faces);
    m_Normals.fill(vec3_t(0,0,0), num_faces);
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
      m_Triangles[id_cell] = Triangle(Point(a[0], a[1], a[2]), Point(b[0], b[1], b[2]), Point(c[0], c[1], c[2]));
      vec3_t u = b - a;
      vec3_t v = c - a;
      m_Normals[id_cell] = u.cross(v);
      m_Normals[id_cell].normalise();
    }
    m_TriangleTree.rebuild(m_Triangles.begin(), m_Triangles.end());
    m_TriangleTree.accelerate_distance_queries();
  }

  int count = 0;
  int total_count = 0;
  //real var[DIM];
  //dim_t<DIM> dim;
  for (size_t i_patch = 0; i_patch < m_PatchGrid->getNumPatches(); ++i_patch) {

    //CartesianPatch* patch = dynamic_cast<CartesianPatch*>(m_PatchGrid->getPatch(i_patch));
    CartesianPatch* patch = NULL;
    if (patch) {
      size_t i_m = patch->sizeI() - 1;
      size_t j_m = patch->sizeJ() - 1;
      size_t k_m = patch->sizeK() - 1;

      if ( i_m == 0 || j_m == 0 || k_m == 0) {
         cout << endl << endl << "Error! i-> " << i_m << " j-> " << j_m << " k-> " << k_m << endl;
         return;
      }

      recursiveLevelSet(patch, 0, 0, 0, i_m, j_m, k_m);

    }
    else {
      levelSetPerCell(i_patch, count, total_count);
    }

    /*
      ++count;
      ++total_count;
      if (count >= 10000) {
        cout << total_count << " cells and " << i_patch << " of " << m_PatchGrid->getNumPatches() << " patches" << endl;
        count = 0;
      }
      */

  }
}

template <unsigned int DIM, unsigned int IVAR>
void DiscreteLevelSet<DIM,IVAR>::recursiveLevelSet(CartesianPatch* patch, size_t i_o, size_t j_o, size_t k_o, size_t i_m, size_t j_m, size_t k_m)
{
  real var[DIM];
  dim_t<DIM> dim;

  vector<size_t> index;
  index.push_back( patch->index(i_o, j_o, k_o) );
  index.push_back( patch->index(i_m, j_o, k_o) );
  index.push_back( patch->index(i_o, j_m, k_o) );
  index.push_back( patch->index(i_o, j_o, k_m) );
  index.push_back( patch->index(i_m, j_m, k_o) );
  index.push_back( patch->index(i_m, j_o, k_m) );
  index.push_back( patch->index(i_o, j_m, k_m) );
  index.push_back( patch->index(i_m, j_m, k_m) );

  vector<real>   g;
  for (int i = 0; i != index.size(); ++i) {
    vec3_t x = patch->xyzoCell(index[i]);
    g.push_back( computePointLevelSet(x) );
  }
  cout  << "     g0-> "   << g[0]
      << "     g1-> "   << g[1]
      << "     g2-> "   << g[2]
      << "     g3-> "   << g[3]
      << "     g4-> "   << g[4]
      << "     g5-> "   << g[5]
      << "     g6-> "   << g[6]
      << "     g7-> "   << g[7] << endl;
  vec3_t x_ooo = patch->xyzoCell(index[0]);
  vec3_t x_ijk = patch->xyzoCell(index[index.size()-1]);

  real d = m_Tol*distance(x_ooo, x_ijk);
  if ( fabs(g[0]) < d || fabs(g[1]) < d || fabs(g[2]) < d || fabs(g[3]) < d ||
       fabs(g[4]) < d || fabs(g[5]) < d || fabs(g[6]) < d || fabs(g[7]) < d)
  {
    /*
  cout << "io-> " << i_o << " j_o-> " << j_o << " k_o-> " << k_o
      << " im-> " << i_m << " j_m-> " << j_m << " k_m-> " << k_m
      << " d-> "  << d  << endl
      << "     g0-> "   << g[0]
      << "     g1-> "   << g[1]
      << "     g2-> "   << g[2]
      << "     g3-> "   << g[3]
      << "     g4-> "   << g[4]
      << "     g5-> "   << g[5]
      << "     g6-> "   << g[6]
      << "     g7-> "   << g[7] << endl;
      */
    if ( (i_m-i_o) == 1 && (j_m-j_o) == 1 && (k_m-k_o) == 1) {
      for (int i = 0; i != index.size(); ++i) {
        size_t i_cell = index[i];
        //patch->getVar(dim, 0, i_cell, var);
        dynamic_cast<Patch*>(patch)->getVar(dim, 0, i_cell, var);
        var[IVAR] = g[i];
        //patch->setVar(dim, 0, i_cell, var);
        dynamic_cast<Patch*>(patch)->setVar(dim, 0, i_cell, var);
      }
    }  else
    if ( (i_m - i_o) >= max(j_m-j_o, k_m-k_o)) {
      int i_new = (i_m - i_o)/2 + i_o;
      //if (i_new == 1) {
      //  recursiveLevelSet(patch, i_o  , j_o, k_o, i_o+1, j_m, k_m );
      //  recursiveLevelSet(patch, i_o+1, j_o, k_o, i_m  , j_m, k_m );
      //}
      //else {
        recursiveLevelSet(patch, i_o  , j_o, k_o, i_new, j_m, k_m );
        recursiveLevelSet(patch, i_new, j_o, k_o, i_m  , j_m, k_m );
      //}
    } else
    if ( (j_m - j_o) >= max(i_m-i_o, k_m-k_o)) {
      int j_new = (j_m - j_o)/2 + j_o;
      //if (j_new == 1) {
      //  recursiveLevelSet(patch, i_o, j_o  , k_o, i_m, j_o+1, k_m );
      //  recursiveLevelSet(patch, i_o, j_o+1, k_o, i_m, j_m  , k_m );
      //}
      //else {
        recursiveLevelSet(patch, i_o, j_o  , k_o, i_m, j_new, k_m );
        recursiveLevelSet(patch, i_o, j_new, k_o, i_m, j_m, k_m );
      //}
    } else
    if ( (k_m - k_o) >= max(i_m-i_o, j_m-j_o)) {
      int k_new = (k_m - k_o)/2 + k_o;
      //if (k_new == 1) {
        //recursiveLevelSet(patch, i_o, j_o, k_o  , i_m, j_m, k_o+1 );
        //recursiveLevelSet(patch, i_o, j_o, k_o+1, i_m, j_m, k_m   );
      //}
      //else {
        recursiveLevelSet(patch, i_o, j_o, k_o  , i_m, j_m, k_new );
        recursiveLevelSet(patch, i_o, j_o, k_new, i_m, j_m, k_m   );
      //}
    }
  }
  else {
    interpolate(patch, i_o, j_o, k_o, i_m, j_m, k_m);
  }
}

template <unsigned int DIM, unsigned int IVAR>
void DiscreteLevelSet<DIM,IVAR>::interpolate(CartesianPatch* patch, size_t i_o, size_t j_o, size_t k_o, size_t i_m, size_t j_m, size_t k_m) {
  cout << endl << endl << "      in interpolate" << endl << endl;
  vec3_t x_ooo = patch->xyzoCell( patch->index(i_o, j_o, k_o) );
  vec3_t x_ijk = patch->xyzoCell( patch->index(i_m, j_m, k_m) );

  vec3_t x_m;
  x_m[0] = x_ijk[0] - x_ooo[0];
  x_m[1] = x_ijk[1] - x_ooo[1];
  x_m[2] = x_ijk[2] - x_ooo[2];

  real g_m = computePointLevelSet(x_m);

  real var[DIM];
  dim_t<DIM> dim;
  for (size_t i = i_o; i <= i_m; ++i) {
    for (size_t j = j_o; j <= j_m; ++j) {
      for (size_t k = k_o; k <= k_m; ++k) {
        patch->getVar(dim, 0, i, j, k, var);
        var[IVAR] = g_m;
        patch->setVar(dim, 0, i, j, k, var);
      }
    }
  }

}

template <unsigned int DIM, unsigned int IVAR>
real DiscreteLevelSet<DIM,IVAR>::distance(vec3_t x1, vec3_t x2) {
  return sqrt(x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]);
}

template <unsigned int DIM, unsigned int IVAR>
real DiscreteLevelSet<DIM,IVAR>::computePointLevelSet(vec3_t x)
{
  real g;
  try {
    Point p(x[0], x[1], x[2]);
    TrianglePointAndPrimitiveId result = m_TriangleTree.closest_point_and_primitive(p);
    Point cp = result.first;
    Triangle* T = result.second;
    int id = (T - m_Triangles.begin());
    vec3_t x_snap = vec3_t(cp[0], cp[1], cp[2]);
    vec3_t v = x - x_snap;
    g = v.abs();
    if (v*m_Normals[id] < 0) {
      g *= -1;
    }
  } catch (...) {
    ERROR("cannot compute distance");
  }
  return g;
}

template <unsigned int DIM, unsigned int IVAR>
void DiscreteLevelSet<DIM,IVAR>::levelSetPerCell(size_t i_patch, int& count, int& total_count)
{
  real var[DIM];
  dim_t<DIM> dim;
  for (size_t i_cell = 0; i_cell < m_PatchGrid->getPatch(i_patch)->variableSize(); ++i_cell) {
    vec3_t x = m_PatchGrid->getPatch(i_patch)->xyzoCell(i_cell);

    m_PatchGrid->getPatch(i_patch)->getVar(dim, 0, i_cell, var);
    var[IVAR] = DiscreteLevelSet::computePointLevelSet(x);
    m_PatchGrid->getPatch(i_patch)->setVar(dim, 0, i_cell, var);

    /*
    ++count;
    ++total_count;
    if (count >= 10000) {
      cout << total_count << " cells and " << i_patch << " of " << m_PatchGrid->getNumPatches() << " patches" << endl;
      count = 0;
    }
    */
  }
}

template <unsigned int IVAR>
struct StoredLevelSet
{
  CUDA_DH static real G(GPU_CartesianPatch& patch, size_t i, size_t j, size_t k, size_t i_field = 0)
  {
    return patch.f(i_field, IVAR, i, j, k);
  }

  CUDA_DH static real G(CartesianPatch& patch, size_t i, size_t j, size_t k, size_t i_field = 0)
  {
    return patch.f(i_field, IVAR, i, j, k);
  }

  CUDA_DH static void updateG(CartesianPatch& patch, size_t i, size_t j, size_t k, real G_value, size_t i_field = 0)
  {
    patch.f(i_field, IVAR, i, j, k) = G_value;
  }

  CUDA_DH static void updateG(GPU_CartesianPatch& patch, size_t i, size_t j, size_t k, real G_value, size_t i_field = 0)
  {
    patch.f(i_field, IVAR, i, j, k) = G_value;
  }
};


#endif // DISCRETELEVELSET_H
