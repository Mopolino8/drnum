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
#ifndef CARTESIANPATCH_H
#define CARTESIANPATCH_H

class CartesianPatch;

//#include <cstddef>
//#include <string.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
#include "patch.h"

#ifdef WITH_VTK
#include <QString>
#include <QVector>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkUnstructuredGrid.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#endif

class CartesianPatch : public Patch
{


#include "cartesianpatch_common.h"

private:

  real m_Lx;
  real m_Ly;
  real m_Lz;

protected: // attributes

  size_t m_NumSeekImin;
  size_t m_NumSeekImax;
  size_t m_NumSeekJmin;
  size_t m_NumSeekJmax;
  size_t m_NumSeekKmin;
  size_t m_NumSeekKmax;

  size_t m_NumProtImin;
  size_t m_NumProtImax;
  size_t m_NumProtJmin;
  size_t m_NumProtJmax;
  size_t m_NumProtKmin;
  size_t m_NumProtKmax;

  size_t m_IDonorZoneFirst;
  size_t m_IDonorZoneAfterlast;
  size_t m_JDonorZoneFirst;
  size_t m_JDonorZoneAfterlast;
  size_t m_KDonorZoneFirst;
  size_t m_KDonorZoneAfterlast;

  size_t m_ICoreFirst;
  size_t m_ICoreAfterlast;
  size_t m_JCoreFirst;
  size_t m_JCoreAfterlast;
  size_t m_KCoreFirst;
  size_t m_KCoreAfterlast;

  real m_LimiterEpsilon;

  real m_xCCMin,      m_yCCMin,      m_zCCMin;
  real m_xCCMax,      m_yCCMax,      m_zCCMax;
  real m_xCCInterMin, m_yCCInterMin, m_zCCInterMin;
  real m_xCCInterMax, m_yCCInterMax, m_zCCInterMax;
  real m_Eps, m_EpsDX, m_EpsDY, m_EpsDZ;

protected: // methods

  void computeDeltas();

  virtual void buildBoundingBox();

public: // methods

  /**
   * Constructor
   * @param num_protectlayers number of protection layers, in which no foreign data access is allowed
   * @param num_overlaplayers number of overlap layers, for which to get data from donor neighbour patches
   */
  CartesianPatch(PatchGrid *patch_grid, size_t num_seeklayers = 2, size_t num_addprotectlayers = 0);


  /**
    * Get neighbour cells for a cell on patch.
    * @param l_cell index of cell for which to search neighbour cells.
    * @param l_cell_neighbours vector with neighbour cell indices.
    */
  virtual void cellNeighbours (const size_t& l_cell,
                               vector<size_t>& i_cell_neighbours);


  /**
    * Get neighbour cells for a cell on patch.
    * @param i_cell i-index of cell for which to search neighbour cells.
    * @param i_cell j-index of cell for which to search neighbour cells.
    * @param i_cell k-index of cell for which to search neighbour cells.
    * @param l_cell_neighbours vector with neighbour cell indices.
    */
  void cellNeighbours (const size_t& i_cell,
                       const size_t& j_cell,
                       const size_t& k_cell,
                       vector<size_t>& l_cell_neighbours);


  /**
    * Read mesh data from file
    * @param s_mesh the stream to read from
    * @return true, if successful
    */
  virtual bool readFromFile(istringstream& iss_input);


  /**
    * Write mesh data to file
    * @param s_mesh the stream to write to
    * @return true, if successful
    */
  //virtual bool writeToFile(ifstream& s_mesh);
  virtual bool writeToFile(ofstream&) {BUG; return false;}


  /**
    * Scale patch relative to origin of parental coordsyst.
    * NOTE: Affects reference position and physical patch size.
    * Virtual: base class Patch holds reference position. Derived class holds phys. patch size.
    * @param scfactor scaling factor.
    */
  virtual void scaleRefParental(real scfactor);


  /**
    * Set individual number of seek layers on all six boundaries of the CartesianPatch. This setting is
    * intended to allow exceptions, like interpolations along boundaries or reduced (1D, 2D) computations.
    * @param num_seekImin number of seek layers on I-min side
    * @param num_seekImax number of seek layers on I-max side
    * @param num_seekJmin number of seek layers on J-min side
    * @param num_seekJmax number of seek layers on J-max side
    * @param num_seekKmin number of seek layers on K-min side
    * @param num_seekKmax number of seek layers on K-max side
    */
  void setSeekExceptions(const size_t& num_seekImin, const size_t& num_seekImax,
                         const size_t& num_seekJmin, const size_t& num_seekJmax,
                         const size_t& num_seekKmin, const size_t& num_seekKmax);


  /**
    * Set up the metrics of the block.
    * @param ilength physical size of block in i direction
    * @param jlength physical size of block in j direction
    * @param klength physical size of block in k direction
    */
  void setupMetrics(real ilength, real jlength, real klength);


  /**
    * Set number of cells in i,j,k
    * @param num_i number of cells in i direction
    * @param num_j number of cells in j direction
    * @param num_k number of cells in k direction
    */
  void resize(size_t num_i, size_t num_j, size_t num_k);


  /**
    * Build up regions (seek, protection and core)
    */
  virtual void buildRegions();


  /**
    * Extract set of data seeking cells on the boundary faces of the patch.
    */
  virtual void extractSeekCells();


  /**
     * Compute dependencies "from" a neighbour.
     * receiving patch: "this"
     * donor patch:     neighbour_patch
     * @param i_neighbour index of neighbour patch from which to receive data
     * @return bool indicating any dependency was found
     */
  virtual bool computeDependencies(const size_t& i_neighbour);


  /**
    * Check overlap with a box defined in xyzo.
    * @param box_xyzo_min lower coords of box
    * @param box_xyzo_min upper coords of box
    * @param only_core indicates to analise only core region of patch
    * @return true, if overlap exists
    */
  virtual bool checkBoxOverlap(const vec3_t& box_xyzo_min, const vec3_t& box_xyzo_max,
                               const bool& only_core = true);


  /// @todo shift all addressing stuff into StructuredPatch and inherite CartesianPatch from it
  /// @todo check, if the int - size_t conversion may cause performance issues

  /**
   * @brief Get the indicees (i, j, k) from field index/
   * @param l the field index
   * @param i first Cartesian index  (return reference)
   * @param j second Cartesian index (return reference)
   * @param k third Cartesian index  (return reference)
   * @return the index in the one dimensional data field
   */
  inline void ijk(const size_t& l,
                  size_t& i, size_t& j, size_t& k) const {
    //    div_t divresult;
    //    divresult = div (l, m_NumJ*m_NumK);  // causes compile error (ambiguous)
    //    i = divresult.quot;
    //    divresult = div (divresult.rem, m_NumK);
    //    j = divresult.quot;
    //    k = divresult.rem;
    size_t rest;
    i = l / (m_NumJ*m_NumK);
    rest = l - i*(m_NumJ*m_NumK);
    j = rest / m_NumK;
    k = rest - j*m_NumK;
  }

  /// @todo too many xyzCell methods. Clean up here.
  inline vec3_t xyzCell(const size_t& l) const
  {
    size_t i, j, k;
    ijk(l,
        i, j, k);
    return xyzCell(i, j, k);
  }


  inline vec3_t xyzCell(const size_t& i, const size_t& j, const size_t& k) const
  {
    vec3_t xyz;
    xyzCell(i, j, k,
            xyz[0], xyz[1], xyz[2]);
    return xyz;
  }


  virtual void xyzCell(const size_t& l_cell,
                       real& x_cell, real& y_cell, real& z_cell)
  {
    size_t i, j, k;
    ijk(l_cell,
        i, j, k);
    xyzCell(i, j, k,
            x_cell, y_cell, z_cell);
  }


  inline void xyzCell(const size_t& i, const size_t& j, const size_t& k,
                      real& x, real& y, real& z) const
  {
    x = m_xCCMin + i*m_Dx;
    y = m_yCCMin + j*m_Dy;
    z = m_zCCMin + k*m_Dz;
  }


  virtual void xyzSubCellRaster(const size_t& l_cell, const size_t& lin_mult_res,
                                vector<vec3_t>& subcells,
                                vec3_t& ref_dxyz)
  {
    subcells.clear();

    // i, j, k
    size_t i, j, k;
    ijk(l_cell,
        i, j, k);

    // increment in subcell resolution
    real resolve_factor = 1./float(lin_mult_res);
    real dx_sub = m_Dx * resolve_factor;
    real dy_sub = m_Dy * resolve_factor;
    real dz_sub = m_Dz * resolve_factor;

    // lower corner point in all coords
    real x_low = i*m_Dx + 0.5 * dx_sub;
    real y_low = j*m_Dy + 0.5 * dy_sub;
    real z_low = k*m_Dz + 0.5 * dz_sub;

    // fill subcell coord triples in list
    vec3_t xyz_h;
    for (size_t i_sub = 0; i_sub < lin_mult_res; i_sub ++) {
      xyz_h[0] = x_low + i_sub * dx_sub;
      for (size_t j_sub = 0; j_sub < lin_mult_res; j_sub ++) {
        xyz_h[1] = y_low + j_sub * dy_sub;
        for (size_t k_sub = 0; k_sub < lin_mult_res; k_sub ++) {
          xyz_h[2] = z_low + k_sub * dz_sub;
          subcells.push_back(xyz_h);
        }
      }
    }
    ref_dxyz[0] = m_Dx;
    ref_dxyz[1] = m_Dy;
    ref_dxyz[2] = m_Dz;
  }


  inline int boundingNVecDirection(const size_t& lc) const
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


  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param error true, if (i, j, k) are out of bounds
   * @return the index in the one dimensional data field
   */
  size_t save_index(size_t i, size_t j, size_t k,
                    bool& error) {
    error = false;
    // NOTE: attempt to enter a negative i,j, of k will be seen as out of range,
    //       as size_t(negative) is huge.
    if(i > m_NumI-1) error = true;
    if(j > m_NumJ-1) error = true;
    if(k > m_NumK-1) error = true;
    return i*m_NumJ*m_NumK + j*m_NumK + k;
  }

//  size_t save_index(int i, int j, int k,
//                    bool& error) {
//    size_t si = i;  // avoid vast compiler warnings
//    size_t sj = j;
//    size_t sk = k;
//    error = false;
//    if(i < 0) error = true;
//    if(j < 0) error = true;
//    if(k < 0) error = true;
//    if(si > m_NumI-1) error = true;
//    if(sj > m_NumJ-1) error = true;
//    if(sk > m_NumK-1) error = true;
//    return si*m_NumJ*m_NumK + sj*m_NumK + sk;
//  }

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


  // new since Sun Jun 30 04:15:41 CEST 2013
  virtual bool computeCCDataInterpolCoeffs_V1 (real x, real y, real z,
                                               WeightedSet<real>& w_set);

  // new since Sun Jun 30 04:15:41 CEST 2013
  bool linearCCInterpolCoeffs (const real& s_test,
                               const size_t& ijk_donor_first, const size_t& ijk_donor_afterlast,
                               const real& ds, const size_t& num,
                               size_t& cell_ref, real& s_coeff, int& sector);

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

  /** Compute ...
    * @return smallest characteristic length.
    */
  virtual real computeMinChLength();

  /**
   * Write patch data to an individual files for this patch only.
   * Example: Let m_myindex be 777 and write to data file with base name "calc/mydata" at counter level 55.
   * Files to write to: calc/mydata_ip000777_000055.vtr
   * @param base_data_filename base data filename relative to cwd.
   * @param count discrete counter (usually time counter).
   */
  virtual void writeData(QString base_data_filename, int count);


  /**
   * Collect neighbour cells of a cell. Version for neighbour cells, that share a
   * common face with cell.
   * @param l_cell index of the cell for which to search neighbour cells.
   * @param l_cell_neighbours vector with neighbour cells indices (ret. ref.)
   */
  virtual void cellOverFaceNeighbours(const size_t& l_c,
                                      vector<size_t>& l_cell_neighbours);


  /**
   * Collect neighbour cells of a cell. Version for neighbour cells, that share at
   * least one common node with cell.
   * @param l_cell index of the cell for which to search neighbour cells.
   * @param l_cell_neighbours vector of neighbour cells of l_cell (ret. ref.)
   */
  virtual void cellOverNodeNeighbours(const size_t& l_c,
                                      vector<size_t>& l_cell_neighbours);


  /**
   * Collect neighbour cells of a cell. Version for neighbour cells, that share a
   * common face with cell.
   * @param l_cell index of the cell for which to search neighbour cells.
   * @param l_cell_neighbours vector with neighbour cells indices (ret. ref.)
   */
  void cellOverFaceNeighbours(const size_t& i_c,
                              const size_t& j_c,
                              const size_t& k_c,
                              vector<size_t>& l_cell_neighbours);


  /**
   * Collect neighbour cells of a cell. Version for neighbour cells, that share at
   * least one common node with cell.
   * @param l_cell index of the cell for which to search neighbour cells.
   * @param l_cell_neighbours vector of neighbour cells of l_cell (ret. ref.)
   */
  void cellOverNodeNeighbours(const size_t& i_c,
                              const size_t& j_c,
                              const size_t& k_c,
                              vector<size_t>& l_cell_neighbours);


#ifdef WITH_VTK
  virtual vtkSmartPointer<vtkDataSet> createVtkDataSet(size_t i_field, const PostProcessingVariables& proc_vars);
  virtual vtkSmartPointer<vtkUnstructuredGrid> createVtkGridForCells(const list<size_t> &cells);
#endif

};

#endif // CARTESIANPATCH_H
