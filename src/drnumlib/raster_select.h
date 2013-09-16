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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef RASTER_H
#define RASTER_H

#include <cstddef>
#include <vector>
#include "blockcfd.h"
#include "weightedset.h"
#include "math/coordtransformvv.h"
//#include "List.hh"

//class Raster;

//#ifdef WITH_VTK
//#include <QString>
//#include <QVector>
//#include <vtkSmartPointer.h>
//#include <vtkRectilinearGrid.h>
//#include <vtkXMLRectilinearGridWriter.h>
//#include <vtkFloatArray.h>
//#include <vtkCellData.h>
//#endif

class Raster : public Raster
{

protected: // attributes

//sh   size_t m_NumI;
//sh   size_t m_NumJ;
//sh   size_t m_NumK;

  real m_Xo;  ///< xo-pos of reference point in intertial coords
  real m_Yo;  ///< yo-pos of reference point in intertial coords
  real m_Zo;  ///< yo-pos of reference point in intertial coords

//cart   real m_Uxo;
//cart   real m_Uyo;
//cart   real m_Uzo;

//cart   real m_Vxo;
//cart   real m_Vyo;
//cart   real m_Vzo;

//cart   real m_Wxo;
//cart   real m_Wyo;
//cart   real m_Wzo;

  real m_Lx;   ///< size in x-direction (own coord syst)
  real m_Ly;   ///< size in y-direction (own coord syst)
  real m_Lz;   ///< size in z-direction (own coord syst)

//cart  real m_Dx;
//cart  real m_Dy;
//cart  real m_Dz;
//cart  real m_InvDx;
//cart  real m_InvDy;
//cart  real m_InvDz;

  /// @todo change to vec3_t later
  //  SVec3 m_xyzCCMin;      ///< Cell center coords of lowest address cell
  //  SVec3 m_xyzInterMin;   ///< lower CC coord-limits to access data from foreign patches
  //  SVec3 m_xyzInterMax;   ///< upper CC coord-limits to access data from foreign patches

  /** @todo attention: m_NumProtectLayers has been defined in "patch". If ever "patch" will
    * be derived from this raster pattern, the following attributes will cause an inheritance conflict:
    * m_transformInertial2This, m_bbox_xyzo_min, m_bbox_xyzo_max, m_NumProtectLayers. */
  CoordTransformVV m_transformInertial2This; ///< transformation matrix to transform intertial coords into system of "this"
  vec3_t m_bbox_xyzo_min;                    ///< lowest coordinates of smallest box around patch in inertial coords.
  vec3_t m_bbox_xyzo_max;                    ///< highest coordinates of smallest box around patch in inertial coords.
  size_t m_NumProtectLayers;                 ///< Number of protection layers

//cart  real m_xCCMin,      m_yCCMin,      m_zCCMin;
//cart  real m_xCCMax,      m_yCCMax,      m_zCCMax;
//cart  real m_xCCInterMin, m_yCCInterMin, m_zCCInterMin;
//cart  real m_xCCInterMax, m_yCCInterMax, m_zCCInterMax;
  real m_Eps, m_EpsDX, m_EpsDY, m_EpsDZ;  ///< precision thresholds
  /** @todo We will need prec limits thoughout the program. It would be usefull to provide basic values depending on data
    *       type used, e.g. float basicEps = 1.e-6 and double basicEps = 1.e-12   This allows to set thresholds independent from
    *       actuial real type in use. */

 // List* m_CellLink;  ///< Link point for FVMOUSE::List derived classes

protected: // methods

//cart  void computeDeltas();


public: // methods

  Raster();
//cart  void setupAligned(real xo1, real yo1, real zo1, real xo2, real yo2, real zo2);
//sh  void resize(size_t num_i, size_t num_j, size_t num_k);
  void setNumProtectLayers(size_t num_protectlayers);  ///< Set number of protection layers

//sh  size_t sizeI() { return m_NumI; }
//sh  size_t sizeJ() { return m_NumJ; }
//sh  size_t sizeK() { return m_NumK; }
  size_t sizeL() { return m_NumI * m_NumJ * m_NumK; }

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @return the index in the one dimensional data field
   */
//sh  size_t index(int i, int j, int k) { return i*m_NumJ*m_NumK + j*m_NumK + k; }

  /**
   * @brief Get the indicees (i, j, k) from field index/
   * @param l the field index
   * @param i first  index  (return reference)
   * @param j second  index (return reference)
   * @param k third  index  (return reference)
   * @return the index in the one dimensional data field
   */
//sh  inline void ijk(const size_t& l,
//sh                  size_t& i, size_t& j, size_t& k) const;

  /**
   * @brief Get cell center from 1D field index/
   * @param l the field index
   * @return the vector xyz of cell center coordinates
   */
//cart   inline vec3_t xyzCell(const size_t& l) const;

  /**
   * @brief Get cell center from indicees (i,j,k)/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @return the vector xyz of cell center coordinates
   */
//cart   inline vec3_t xyzCell(const size_t& i, const size_t& j, const size_t& k) const;

  /**
   * @brief Get cell center from indicees (i,j,k)/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @param x x-coord of cell center (return reference)
   * @param y y-coord of cell center (return reference)
   * @param z z-coord of cell center (return reference)
   */
//cart   inline void xyzCell(const size_t& i, const size_t& j, const size_t& k,
//cart                       real& x, real& y, real& z) const;

  /**
   * Build a bounding box in inertial coordinates around the  raster.
   * NOTE: Raster coords may be non-aligned to the xyzo-system.
   */
  virtual void buildBoundingBox();

  /**
   * @brief Get normal direction of nearest boundary of the raster
   * @param lc the 1D cell field index
   * @return an integer 0, 1 or 2 for the directions i,j,k respectively
   */
//cart   inline int boundingNVecDirection(const size_t& lc) const;

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @param error true, if (i, j, k) are out of bounds
   * @return the index in the one dimensional data field
   */
//sh   size_t save_index(int i, int j, int k,
//sh                     bool& error);

  /**
   * Check if an (i,j,k) triple is inside the patch.
   * Attention only the upper limit will be checked (unsigned data type).
   * @param i first index
   * @param j second index
   * @param k third index
   * @return true if it is a valid (i,j,k) triple, false otherwise
   */
//sh   bool checkRange(size_t i, size_t j, size_t k);


//cart   real dx() { return m_Dx; }
//cart   real dy() { return m_Dy; }
//cart   real dz() { return m_Dz; }
//cart   real dV() { return m_Dx*m_Dy*m_Dz; }
//cart   real idx() { return m_InvDx; }
//cart   real idy() { return m_InvDy; }
//cart   real idz() { return m_InvDz; }

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
//cart   void computeCCInterpolWeights(real& x, real& y, real& z,
//cart                                 size_t& ic_ref, size_t& jc_ref, size_t& kc_ref,
//cart                                 real* w_octet);

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
//cart  void computeNCInterpolWeights(real& x, real& y, real& z,
//cart                                size_t& in_ref, size_t& jn_ref, size_t& kn_ref,
//cart                                real* w_octet);

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
//cart  void computeInterpolWeights(real& x, real& y, real& z,
//cart                              real& x_ref, real& y_ref, real& z_ref,
//cart                              real* w_octet);

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
//cart  bool xyzToRefInterCell(real x, real y, real z,
//cart			 size_t& i_ref, size_t& j_ref, size_t& k_ref);

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
//cart  bool xyzToRefCell(real x, real y, real z,
//cart                    size_t& ic_ref, size_t& jc_ref, size_t& kc_ref);

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
//cart  bool xyzToRefNode(real x, real y, real z,
//cart                    size_t& in_ref, size_t& jn_ref, size_t& kn_ref);

};

//sh inline bool Raster::checkRange(size_t i, size_t j, size_t k)
//sh {
//sh   if (i >= sizeI() || j >= sizeJ() || k >= sizeK()) {
//sh     return false;
//sh   }
//sh   return true;
//sh }


//sh inline void Raster::ijk(const size_t& l,
//sh                                  size_t& i, size_t& j, size_t& k) const {
//sh   size_t rest;
//sh   i = l / (m_NumJ*m_NumK);
//sh   rest = l - i*(m_NumJ*m_NumK);
//sh   j = rest / m_NumK;
//sh   k = rest - j*m_NumK;
//sh }


//cart inline vec3_t Raster::xyzCell(const size_t& l) const
//cart {
//cart   size_t i, j, k;
//cart   ijk(l,
//cart       i, j, k);
//cart   return xyzCell(i, j, k);
//cart }


//cart inline vec3_t Raster::xyzCell(const size_t& i, const size_t& j, const size_t& k) const
//cart {
//cart   vec3_t xyz;
//cart   xyzCell(i, j, k,
//cart           xyz[0], xyz[1], xyz[2]);
//cart   return xyz;
//cart }


//cart inline void Raster::xyzCell(const size_t& i, const size_t& j, const size_t& k,
//cart                                      real& x, real& y, real& z) const
//cart {
//cart   x = m_xCCMin + i*m_Dx;
//cart   y = m_yCCMin + j*m_Dy;
//cart   z = m_zCCMin + k*m_Dz;
//cart }


//cart inline int Raster::boundingNVecDirection(const size_t& lc) const
//cart {
//cart   int direction;
//cart   size_t i, j, k;
//cart   size_t min_index_dist, min_index_dist_h;
//cart   ijk(lc,
//cart       i, j, k);
//cart   // prime with i-direction: 0
//cart   min_index_dist = min(i, (m_NumI-1-i));
//cart   direction = 0;
//cart   // check j-direction: 1
//cart   min_index_dist_h = min(j, (m_NumJ-1-j));
//cart   if(min_index_dist > min_index_dist_h) {
//cart     min_index_dist = min_index_dist_h;
//cart     direction = 1;
//cart   }
//cart   min_index_dist_h = min(k, (m_NumK-1-k));
//cart   if(min_index_dist > min_index_dist_h) {
//cart     min_index_dist = min_index_dist_h;
//cart     direction = 2;
//cart   }
//cart   return direction;
//cart }


//sh inline size_t Raster::save_index(int i, int j, int k,
//sh                                           bool& error)
//sh {
//sh   size_t si = i;  // avoid vast compiler warnings
//sh   size_t sj = j;
//sh   size_t sk = k;
//sh   error = false;
//sh   if(i < 0) error = true;
//sh   if(j < 0) error = true;
//sh   if(k < 0) error = true;
//sh   if(si > m_NumI-1) error = true;
//sh   if(sj > m_NumJ-1) error = true;
//sh   if(sk > m_NumK-1) error = true;
//sh   return si*m_NumJ*m_NumK + sj*m_NumK + sk;
//sh }


#endif // RASTER_H
