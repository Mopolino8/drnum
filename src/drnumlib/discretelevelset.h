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

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>

#include <vtkPolyData.h>

#include "patchgrid.h"
#include "postprocessingvariables.h"

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

  PatchGrid*           m_PatchGrid;
  unsigned int         m_VarIndex;
  //SegmentTree          m_SegmentTree;
  //vtkPolyData*         m_PolyData;
  //QVector<vtkIdType>   m_Tri2Grid;
  //QVector<double>      m_Radius;      ///< Surface radius for mesh resolution.


public:

  DiscreteLevelSet(PatchGrid *patch_grid, int var_index, QString stl_file_name);


public: // types

  template <unsigned int IVAR>
  class PlotVars : public PostProcessingVariables
  {
  public:

    virtual int numScalars() const { return 1; }
    virtual int numVectors() const { return 0; }

    virtual string getScalarName(int) const { return "G"; }
    virtual string getVectorName(int) const { BUG; }
    virtual real   getScalar(int, real* var) const { return var[IVAR]; }
    virtual vec3_t getVector(int, real*)     const { BUG; return vec3_t(0,0,0); }
  };

};

#endif // DISCRETELEVELSET_H
