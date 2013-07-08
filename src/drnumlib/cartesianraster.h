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
#ifndef CARTESIANRASTER_H
#define CARTESIANRASTER_H

//#include <cstddef>
#include "structuredhexraster.h"

//class CartesianRaster;

//#ifdef WITH_VTK
//#include <QString>
//#include <QVector>
//#include <vtkSmartPointer.h>
//#include <vtkRectilinearGrid.h>
//#include <vtkXMLRectilinearGridWriter.h>
//#include <vtkFloatArray.h>
//#include <vtkCellData.h>
//#endif

class CartesianRaster : public StructuredHexRaster
{

protected: // attributes

  real m_Uxo;
  real m_Uyo;
  real m_Uzo;

  real m_Vxo;
  real m_Vyo;
  real m_Vzo;

  real m_Wxo;
  real m_Wyo;
  real m_Wzo;


  real m_Dx;
  real m_Dy;
  real m_Dz;
  real m_InvDx;
  real m_InvDy;
  real m_InvDz;

  real m_xCCMin,      m_yCCMin,      m_zCCMin;
  real m_xCCMax,      m_yCCMax,      m_zCCMax;
  real m_xCCInterMin, m_yCCInterMin, m_zCCInterMin;
  real m_xCCInterMax, m_yCCInterMax, m_zCCInterMax;

protected: // methods

  void computeDeltas();


public: // methods

  CartesianRaster();

  virtual void resize(size_t num_i, size_t num_j, size_t num_k);

  void setupAligned(real xo1, real yo1, real zo1, real xo2, real yo2, real zo2);


  /**
   * @brief Get cell center from 1D field index/
   * @param l the field index
   * @return the vector xyz of cell center coordinates
   */
  inline vec3_t xyzCell(const size_t& l) const;

  /**
   * @brief Get cell center from indicees (i,j,k)/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @return the vector xyz of cell center coordinates
   */
  inline vec3_t xyzCell(const size_t& i, const size_t& j, const size_t& k) const;

  /**
   * @brief Get cell center from indicees (i,j,k)/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @param x x-coord of cell center (return reference)
   * @param y y-coord of cell center (return reference)
   * @param z z-coord of cell center (return reference)
   */
  inline void xyzCell(const size_t& i, const size_t& j, const size_t& k,
                      real& x, real& y, real& z) const;

  /**
   * Build a bounding box in inertial coordinates around the  CartesianRaster.
   * NOTE: CartesianRaster coords may be non-aligned to the xyzo-system.
   */
  virtual void buildBoundingBox();

  /**
   * @brief Get normal direction of nearest boundary of the raster
   * @param lc the 1D cell field index
   * @return an integer 0, 1 or 2 for the directions i,j,k respectively
   */
  inline int boundingNVecDirection(const size_t& lc) const;


  real dx() { return m_Dx; }
  real dy() { return m_Dy; }
  real dz() { return m_Dz; }
  real dV() { return m_Dx*m_Dy*m_Dz; }
  real idx() { return m_InvDx; }
  real idy() { return m_InvDy; }
  real idz() { return m_InvDz; }

  /**
   * Set up interpolation methods for giving data to foreign patches.
   * Example: Build up Split- or Octrees for search operations, etc ... depending on patch type.
   * @param num_protection number of overlap cell layers in which no data access is permissible.
   */
  virtual void setupInterpolators();

  /**
   * Get data interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCDataInterpolCoeffs(real x, real y, real z,
                                           WeightedSet<real>& w_set);

  /**
   * Get directional derivative (grad*n) interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param nx the x-component of directional vector in the coords of the present patch
   * @param ny the y-component of directional vector in the coords of the present patch
   * @param nz the z-component of directional vector in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCGrad1NInterpolCoeffs(real x, real y, real z,
					     real nx, real ny, real nz,
                                             WeightedSet<real>& w_set);

  /**
   * Compute a set of CC interpolation weights for a given point (x,y,z) in inverse mesh.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param ic_ref i-index of reference cell
   * @param jc_ref j-index of reference cell
   * @param kc_ref k-index of reference cell
   * @param w_octet eight values, that represent interpolation weights for cells in intervall (return reference)
   */
  void computeCCInterpolWeights(real& x, real& y, real& z,
                                size_t& ic_ref, size_t& jc_ref, size_t& kc_ref,
                                real* w_octet);

  /**
   * Compute a set of NC interpolation weights for a given point (x,y,z) in natural mesh.
   * NOTE: mesh cell formed by nodes in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param in_ref i-index of reference node
   * @param jn_ref j-index of reference node
   * @param kn_ref k-index of reference node
   * @param w_octet eight values, that represent interpolation weights for cells in intervall (return reference)
   */
  void computeNCInterpolWeights(real& x, real& y, real& z,
                                size_t& in_ref, size_t& jn_ref, size_t& kn_ref,
                                real* w_octet);

  /**
   * Compute a set of interpolation weights for a given point (x,y,z) in a box.
   * NOTE: box corner coords in 3D intervall
   * ([x_ref..x_ref+m_Dx], [y_ref..y_ref+m_Dy], [z_ref..z_ref+m_Dz]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param in_ref i-index of reference node
   * @param jn_ref j-index of reference node
   * @param kn_ref k-index of reference node
   * @param w_octet eight values, that represent interpolation weights for cells in intervall (return reference)
   */
  void computeInterpolWeights(real& x, real& y, real& z,
                              real& x_ref, real& y_ref, real& z_ref,
                              real* w_octet);

  /**
   * Find inverse mesh cell, in which the point (x,y,z) resides.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param i_ref i-index of reference cell (return reference)
   * @param j_ref j-index of reference cell (return reference)
   * @param k_ref k-index of reference cell (return reference)
   * @return boolean, indicating point (x,y,z) is inside core region
   */
  bool xyzToRefInterCell(real x, real y, real z,
                         size_t& i_ref, size_t& j_ref, size_t& k_ref);

  /**
   * Find inverse mesh cell, in which the point (x,y,z) resides.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param ic_ref i-index of reference cell (return reference)
   * @param jc_ref j-index of reference cell (return reference)
   * @param kc_ref k-index of reference cell (return reference)
   * @return boolean, indicating point (x,y,z) is inside core region
   */
  bool xyzToRefCell(real x, real y, real z,
                    size_t& ic_ref, size_t& jc_ref, size_t& kc_ref);

  /**
   * Find cell, in which the point (x,y,z) resides.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param in_ref i-index of reference node (return reference)
   * @param jn_ref j-index of reference node (return reference)
   * @param kn_ref k-index of reference node (return reference)
   * @return boolean, indicating point (x,y,z) is inside core region
   */
  bool xyzToRefNode(real x, real y, real z,
                    size_t& in_ref, size_t& jn_ref, size_t& kn_ref);

};

inline vec3_t CartesianRaster::xyzCell(const size_t& l) const
{
  size_t i, j, k;
  ijk(l,
      i, j, k);
  return xyzCell(i, j, k);
}


inline vec3_t CartesianRaster::xyzCell(const size_t& i, const size_t& j, const size_t& k) const
{
  vec3_t xyz;
  xyzCell(i, j, k,
          xyz[0], xyz[1], xyz[2]);
  return xyz;
}


inline void CartesianRaster::xyzCell(const size_t& i, const size_t& j, const size_t& k,
                                     real& x, real& y, real& z) const
{
  x = m_xCCMin + i*m_Dx;
  y = m_yCCMin + j*m_Dy;
  z = m_zCCMin + k*m_Dz;
}


inline int CartesianRaster::boundingNVecDirection(const size_t& lc) const
{
  int direction;
  size_t i, j, k;
  size_t min_index_dist, min_index_dist_h;
  ijk(lc,
      i, j, k);
  // prime with i-direction: 0
  min_index_dist = min(i, (m_NumI-1-i));
  direction = 0;
  // check j-direction: 1
  min_index_dist_h = min(j, (m_NumJ-1-j));
  if(min_index_dist > min_index_dist_h) {
    min_index_dist = min_index_dist_h;
    direction = 1;
  }
  min_index_dist_h = min(k, (m_NumK-1-k));
  if(min_index_dist > min_index_dist_h) {
    min_index_dist = min_index_dist_h;
    direction = 2;
  }
  return direction;
}


#endif // CARTESIANRASTER_H
