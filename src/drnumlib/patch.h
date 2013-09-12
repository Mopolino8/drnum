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
#ifndef PATCH_H
#define PATCH_H

#include "blockcfd.h"

#include <iostream>
#include <sstream>
#include <string>
#include <list>

#include <QString>

#include <vtkDataSet.h>


#include "weightedset.h"

struct donor_t;
class Patch;
class PatchGrid;
class GPU_Patch;

#include "intercoeffpad.h"
#include "intercoeffws.h"
#include "math/coordtransformvv.h"
#include "codestring.h"
#include "postprocessingvariables.h"
#include "donor_t.h"

#ifdef CUDA
#include "cudatools.h"
#endif

#ifdef WITH_VTK
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#endif

/** @todo
 *    proposed coord naming convention:
 *        (Xo,Yo,Zo), (xo,yo,zo) or vec3_t XYZo, xyzo     : coords in syst. of origin
 *        (X,Y,Z)   , (x,y,z)    or vec3_t XYZ , xyz      : coords in syst. of "this" patch
 *        (XX,YY,ZZ), (xx,yy,zz) or vec3_t XXYYZZ, xxyyzz : coords in syst. of any foreign patch
 *
 */

/**
 * @todo Build new heritage tree for patches. Proposed as follows:
 *  base:        Patch
 *  addr. level: StructuredHexPatch, StructuredPrisPatch, SemistructPrisPatch, SemistructHexPatch, UnstructuredPatch
 *  geo. level:  StructuredHexPatch  => CartesianPatch, CartesianStretchedPatch, GeneralStructuredHexPatch, etc...
 *               StructuredPrisPatch => GeneralStructuredPrisPatch, ... (and geometrically simplified types)
 *               SemistructPrisPatch => ... (unstructured 2D-hex grid, structured in 3rd dimension)
 *               SemistructHexPatch  => ... (unstructured 2D-hex grid, structured in 3rd dimension)
 *               UnstructuredPatch   => ... (probably no variants, hybrid cells)
 *
 *  Remark: at present I think, we should NOT differentiate patch alignments to (xo,yo,zo) from others. Might be
 *          smarter to define coalignments of neighbouring blocks.
 *
 */
//=======
#include "transformation.h"    /// @todo merge: kept for compatibility
//>>>>>>> master

class Patch
{

  friend class GPU_Patch;

#include "patch_common.h"

  void  allocateData();


private: // attributes

  real*    m_GpuData;
  size_t*  m_GpuReceivingCellIndicesConcat;
  size_t*  m_GpuReceivingCellIndicesUnique;
  donor_t* m_GpuDonors;
  size_t*  m_GpuDonorIndexConcat;
  real*    m_GpuDonorWeightConcat;
  bool     m_GpuDataSet;


protected: // attributes

  size_t m_MyIndex;           ///< Index of patch in sequence of PatchGrid::m_patches. Optional setting.
  size_t m_mypatchtype;       ///< Type-code as defined by derived patches. Example 1001: CartesianPatch, etc...
  string m_patchcomment;      ///< Optional comment for this patch to be read/written in patchgrid files
  CodeString m_solvercodes;   ///< String containing infos for choice of field-fluxes, RB-fluxes, sources, ...

  // reference position
  /// @todo redundant info: ommit m_Xo,m_Yo,m_Zo as it is contained in m_transformInertial2This
  real m_Xo;                  ///< X-coord of reference position in parental coordinates
  real m_Yo;                  ///< Y-coord of reference position in parental coordinates
  real m_Zo;                  ///< Z-coord of reference position in parental coordinates
  // settings
  bool m_InterpolateData;     ///< Flag indicates wether to interpolate data on interpatch transfers
  //  bool m_InterpolateGrad1N;   ///< Flag indicates wether to interpolate directed gradients on interpatch transfers
  bool m_TransferPadded;      ///< Flag indicates wether to transfer donor data in padded versions with "InterCoeffPad".
  bool m_SeekExceptions ;     ///< Flag indicates a seek exception: the number of seek layers is not constant.
  size_t m_NumSeekLayers;     ///< number of boundary cell layers, for which to get data from donor neighbour patches
  size_t m_NumAddProtectLayers; ///< additional number of boundary protected layers, in which no interpol access from other patches is allowed

  // IO-scaling, position, orientation
  real m_ioscale;                            ///< general scale factor of mesh related to io-values
  CoordTransformVV m_TransformInertial2This; ///< transformation matrix to transform intertial coords into system of "this"

  // bounding box
  vec3_t m_BBoxXYZoMin;   ///< lowest coordinates of smallest box around patch in inertial coords.
  vec3_t m_BBoxXYZoMax;   ///< highest coordinates of smallest box around patch in inertial coords.
  bool   m_BBoxOk;        ///< flag indicating wether the bounding box is available

  // intermediate variables
  bool m_receiveCells_OK; ///< Flag indicating that receive_cells have been extracted yet. Might be set false on mesh changes.

  // lists related to receiving cells in overlap layers
  vector<size_t> m_ReceiveCells;           ///< cells of "this", expecting to get data from any donor neighbour
  vector<size_t> m_receive_cell_data_hits;  ///< number of contributing patches for data, note indexing as receive_cells
  //  vector<size_t> m_receive_cell_grad1N_hits;///< number of contributing patches for grad1N, note indexing as receive_cells

  // lists related to neighbouring donor patches
  vector<pair<Patch*, CoordTransformVV> > m_neighbours; ///< neighbouring donor patches and coord transformation. NOTE: receiving patches not stored.
  vector<InterCoeffPad> m_InterCoeffData;               ///< Interpolation coefficient lists for data
  //  vector<InterCoeffPad> m_InterCoeffGrad1N;             ///< Interpolation coefficient lists for 1st directed gradients
  vector<InterCoeffWS> m_InterCoeffData_WS;             ///< same as m_InterCoeffData, but with WeightedSets (CPU only)
  //  vector<InterCoeffWS> m_InterCoeffGrad1N_WS;           ///< same as m_InterCoeffGrad1N, but with WeightedSets (CPU only)
  // postponed  vector<InterCoeff*> m_InterCoeffGrad2N;    ///< Interpolation coefficient lists for 2nd directed gradients

  vector<size_t> m_VectorVarIndices; ///< array containing the starting indices of all vectorial variables (e.g. velocity, momentum, ...)


protected: // methods

  // internal data handling
  void  deleteData();
  void  resize(size_t variable_size);

  void setTransformation(Transformation t) { m_Transformation = t; }  /// @todo keep for compatibility, prefer CoordTransformVV later

  // geometry
  virtual void buildBoundingBox()=0;

  // external data handling
  //virtual void extractReceiveCells() { BUG; } ///< Extract cells in overlap layers (m_receive_cells)
  // NEW_SEEK_EXCEPTION
  //virtual void extractSeekCells() { BUG; } ///< Extract seeking cells in overlap layers (m_receive_cells)

  /**
    * Build up regions (seek, protection and core)
    */
  virtual void buildRegions() {BUG;}


public: // methods


  /**
   * Constructor
   * @param patch_grid the grid this patch belonges to
   * @param num_seeklayers default number of seeking element layers
   * @param num_addprotectlayers default number of additional protection layers
   */
  Patch(PatchGrid* patch_grid, size_t num_seeklayers = 2, size_t num_addprotectlayers = 0);


  /**
    * Build direct transfer lists for donor->this relation upon previously build
    * vector<InterCoeffPad> m_InterCoeffData
    */
  void buildDonorTransferData();


  /**
    * Get / compute:
    * number of cells in "this" receiving data from neighbour, stride for receiving
    * number of cells in "neighbour" to be served fron "this", stride for serving
    * @param neighbour pointer to neighbour patch to analyse
    * @param vice_exist true, if a donor relation "neighbour"->"this" exists (return reference)
    * @param num_receiving number receiving cells in "this" (return reference)
    * @param receive_stride stride, donor cells per rec. cell (return reference)
    * @param versa_exist true, if a donor relation "this"->"neighbour" exists (return reference)
    * @param num_serving number of cells in neighbour, being served by "this"  (return reference)
    * @param serve_stride stride, donor cells in this per rec. cell in neighbour (return reference)
    * @param vice_versa indicates, wether to compute the "versa" relation relation
    */
  void diagnoseViceVersaDependencies(Patch* neighbour,
                                     bool& vice_exist, size_t& num_receiving, size_t& receive_stride,
                                     bool& versa_exist, size_t& num_seving, size_t& serve_stride,
                                     const bool& vice_versa = true);


  /** @todo decide wether it is more useful to make PatchGrid a friend class
    * and keep extractSeekCells() and compactReceiveCellLists() protected. */
  /**
    * Extract set of data seeking cells on the boundary faces of the patch.
    */
  virtual void extractSeekCells() {BUG;}


  /**
    * Clean up m_receive_cells and hit counter lists
    */
  void compactReceiveCellLists();


  /**
    * Access coordinates of a cell in local xyz-system.
    * @param l_cell index of cell
    * @param x_cell x-coord. of cell center (return reference)
    * @param y_cell y-coord. of cell center (return reference)
    * @param z_cell z-coord. of cell center (return reference)
    */
  //virtual void xyzCell(const size_t& l_cell,
  //                     real& x_cell, real& y_cell, real& z_cell){BUG;};
  virtual void xyzCell(const size_t&,
                       real&, real&, real&){BUG;};



  /**
    * Access coordinates of a cell in inertial xyzo-system.
    * @param l_cell index of cell
    * @param xo_cell xo-coord. of cell center (return reference)
    * @param yo_cell yo-coord. of cell center (return reference)
    * @param zo_cell zo-coord. of cell center (return reference)
    */
  void xyzoCell(const size_t& l_cell,
                real& xo_cell, real& yo_cell, real& zo_cell);


  /**
    * Create a subcell resolution raster of a cell. Coord-syst. of "this".
    * @param l_cell the cell index
    * @param lin_mult_res linear resolutiojn multiplyer. NOTE 3D!
    * @param xxyyzz_subcells vector of subcell coordinate tiples (return reference)
    * @param ref_dxxyyzz reference cell size dxx,dyy,dzz as vec3_t (return reference)
    */
  //  virtual void xxyyzzSubCellRaster(const size_t& l_cell, const size_t& lin_mult_res,
  //                                   vector<vec3_t>& xxyyzz_subcells,
  //                                   vec3_t& ref_dxxyyzz){BUG;};
  virtual void xxyyzzSubCellRaster(const size_t&, const size_t&,
                                   vector<vec3_t>&,
                                   vec3_t&){BUG;};


  /**
    * Create a subcell resolution raster of a cell. Coord-syst. xyzo.
    * @param l_cell the cell index
    * @param lin_mult_res linear resolutiojn multiplyer. NOTE 3D!
    * @param xyzo_subcells vector of subcell coordinate tiples (return reference)
    * @param ref_dxyzo reference cell size x,y,z as vec3_t (return reference)
    */
  void xyzoSubCellRaster(const size_t& l_cell, const size_t& lin_mult_res,
                         vector<vec3_t>& xyzo_subcells,
                         vec3_t& ref_dxyzo);


  /**
    * Data access from all donor patches from direct data lists
    * ATTENTION: DOES NOT TURN ANY VECTORIAL VARIABLES !!
    * @param field the field, for which all variables are transfered
    */
  void accessDonorDataDirect(const size_t &field);


  void  setNumberOfFields(size_t num_fields) { m_NumFields = num_fields; }
  void  setNumberOfVariables(size_t num_variables) { m_NumVariables = num_variables; }


  /**
   * Compute the difference between two variables. This is intended to be used for convergence monitoring.
   * @param i_field1 index of the first field
   * @param i_var1 index of the first variable
   * @param i_field2 index of the first field
   * @param i_var2 index of the second variable
   * @param max_norm will hold the maximal absolute difference
   * @param l2_norm will hold the L2 norm of the difference
   */
  void computeVariableDifference(size_t i_field1, size_t i_var1, size_t i_field2, size_t i_var2, real &max_norm, real &l2_norm);


  /// @todo The following 2 methods might perhaps go, as default is set construcion time (?)
  /**
    * Set default number of additional protection layers between seek and core region
    * @param num_addprotectlayers number of additional protection layers
    */
  void setNumAddProtectLayers(size_t num_addprotectlayers)
  {
    m_NumAddProtectLayers = num_addprotectlayers;
  }


  /**
    * Set default number of element layers, seeking data from other patches.
    * @param num_seeklayers default number of seeking element layers
    */
  void setNumSeekLayers(size_t num_seeklayers)
  {
    m_SeekExceptions = false;
    m_NumSeekLayers = num_seeklayers;
  }


  /**
    * Set interaction with/without data transfers
    * @param interpolate_data bool to cause data interpolation if true
    */
  void setInterpolateData(bool interpolatedata = true)
  {
    m_InterpolateData = interpolatedata;
  }


  /**
    * Set all dependency transfers from any donors to be padded, employing data transfer
    * classes of type "InterCoeffPad".
    * @param trans_padded bool to cause padded data transfers
    */
  void setTransferPadded(bool trans_padded = true)
  {
    m_TransferPadded = trans_padded;
  }


  /**
    * Read mesh data from file
    * @param iss_input the stream to read from
    * @return true, if successful
    */
  virtual bool readFromFile(istringstream& iss_input);


  /**
    * Read solver codes from file
    * @param iss_input the stream to read from
    */
  virtual void readSolverCodes(istringstream& iss_input);


  /**
    * Write mesh data to file
    * @param s_mesh the stream to write to
    * @return true, if successful
    */
  virtual bool writeToFile(ofstream&) { BUG; return true; }
  //  virtual bool writeToFile(ifstream& s_mesh) {return true;};


  /**
   * Write patch data to an individual files for this patch only.
   * Example: Let m_myindex be 777 and write to data file with base name "calc/mydata" at counter level 55.
   * Files to write to: calc/mydata_ip000777_000055.vtr
   * @param base_data_filename base data filename relative to cwd.
   * @param count discrete counter (usually time counter).
   */
  //virtual void writeData(QString base_data_filename, size_t count) {BUG;}
  virtual void writeData(QString, size_t) {BUG;}


  /**
   * @brief Create a vtkDataSet of this Patch.
   * @return a pointer to the vtkDataSet which has been created
   */
#ifdef WITH_VTK
  virtual vtkSmartPointer<vtkDataSet> createVtkDataSet(size_t i_field, const PostProcessingVariables& proc_vars) = 0;
#endif

  /**
   * @brief create a vtkUnstructuredGrid which represents a subset of the whole patch.
   * This method assumes that the concept of cells exists for any kind of patch.
   * @param cells the cells to put into the unstructured grid output
   * @return a pointer to the newly create vtkUnstructuredGrid
   */
#ifdef WITH_VTK
  virtual vtkSmartPointer<vtkUnstructuredGrid> createVtkGridForCells(const list<size_t>& cells) = 0;
#endif



  /**
    * Build up transformation matrix.
    * @param xyzoref position of origin of own coord-syst in parental coord system.
    * @param base_i base vector in i-direction of block, given in parental coords
    * @param base_j base vector in j-direction of block, given in parental coords
    */
  void setupTransformation(vec3_t xyzoref,
                           vec3_t base_i, vec3_t base_j);


  /**
    * Scale patch relative to origin of parental coordsyst.
    * NOTE: Affects reference position and physical patch size.
    * Virtual, base class patch only holds reference position.
    * @param scfactor scaling factor.
    */
  virtual void scaleRefParental(real scfactor)
  {
    m_TransformInertial2This.scaleVector(scfactor);
  }


  /**
    * Set
    @param index the index of the patch in sequence of PatchGrid::m_patches
    */
  void setIndex(size_t patchindex) {m_MyIndex = patchindex;}


  /**
    * Access
    * @return the index of the patch in sequence of PatchGrid::m_patches
    */
  size_t accessIndex() {return m_MyIndex;}


  /**
    * Access
    * @return the number of donor neighbours
    */
  size_t accessNumNeighbours() {return m_neighbours.size();}


  /**
    * Access index of neighbour in sequence of patchgrid::m_Patches
    * @param ii_n the internal counter index of the neighbour in this->m_neighbours
    * @return the pointer to the donor neighbour
    */
  Patch* accessNeighbour(size_t ii_n) {
    return m_neighbours[ii_n].first;
  }


  /**
    * Access index of neighbour in sequence of patchgrid::m_Patches
    * @param ii_n the internal counter index of the neighbour in this->m_neighbours
    * @return the index of the donor neighbour
    */
  size_t accessNeighbourIndex(size_t ii_n) {
    Patch* neighbour = m_neighbours[ii_n].first;
    return neighbour->accessIndex();
  }


  /**
    * Access
    * @return type-code
    */
  size_t accessPatchType() {return m_mypatchtype;}


  /**
    * Set
    * @param patch-comment
    */
  void setPatchComment(string patchcomment) {
    m_patchcomment = patchcomment;
  }


  /**
    * Access
    * @return patch-comment
    */
  string accessPatchComment() {return m_patchcomment;}


  /**
    * Access
    * @return solver-codes
    */
  string accessSolverCodes() {return m_solvercodes;}


  /**
    * Access
    * @return lower coordinates of bounding box
    */
  vec3_t accessBBoxXYZoMin();


  /**
    * Access
    * @return upper coordinates of bounding box
    */
  vec3_t accessBBoxXYZoMax();


  /**
   * Check receiving of any data from a neighbouring donor patch
   * an,d if so, insert it.
   * receiving block: "this"
   * donor block:     neighbour_patch
   * @param neighbour_patch the new donor neighbour patch of "this".
   */
  void insertNeighbour(Patch* neighbour_patch);


  /**
     * Compute dependencies "from" a neighbour.
     * receiving patch: "this"
     * donor patch:     neighbour_patch
     * @param i_neighbour index of neighbour patch from which to receive data
     * @return bool indicating any dependency was found
     */
  virtual bool computeDependencies(const size_t& i_neighbour)=0;


  /**
     * Finalize the computation of dependencies from neighbouring patches.
     * - Reduce contribution weights for receiving cells, being influenced by more than one donor patch.
     * - Transfer data to padded data sets, if required.
     */
  void finalizeDependencies();


  /**
    * Check overlap with a box defined in xyzo.
    * @param box_xyzo_min lower coords of box
    * @param box_xyzo_min upper coords of box
    * @param only_core indicates to analise only core region of patch
    * @return true, if overlap exists
    */
  virtual bool checkBoxOverlap(const vec3_t& box_xyzo_min, const vec3_t& box_xyzo_max,
                               const bool& only_core = true)=0;


  /**
   * Set up interpolation methods for giving data to foreign patches.
   * Example: Build up Split- or Octrees for search operations, etc ... depending on patch type.
   */
  virtual void setupInterpolators()=0;


  /**
   * Get data interpolation coeff-sets.
   * @param xyz coordinates of requesting point in system of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  inline bool computeCCDataInterpolCoeffs(vec3_t xyz,
                                          WeightedSet<real>& w_set) {
    return computeCCDataInterpolCoeffs(xyz[0], xyz[1], xyz[2],
                                       w_set);
  }


  /**
   * Get data interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCDataInterpolCoeffs(real x, real y, real z,
                                           WeightedSet<real>& w_set)=0;

  // new since Sun Jun 30 04:15:41 CEST 2013
  virtual bool computeCCDataInterpolCoeffs_V1 (real x, real y, real z,
                                               WeightedSet<real>& w_set)=0;


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
  inline bool computeCCGrad1NInterpolCoeffs(vec3_t xyz, vec3_t nxyz,
                                            WeightedSet<real>& w_set) {
    return computeCCGrad1NInterpolCoeffs(xyz[0], xyz[1], xyz[2],
                                         nxyz[0], nxyz[1], nxyz[2],
                                         w_set);
  }


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
                                             WeightedSet<real>& w_set)=0;


  /**
    * Data access from all donor patches via m_InterCoeffData_WS
    * NOTE: Non padded version employing WeightedSet pattern.
    * ATTENTION: DOES NOT TURN ANY VECTORIAL VARIABLES !!
    * @param field the field, for which all variables are transfered
    */
  void accessDonorData_WS(const size_t& field);

  /// @todo replace accessDonorData... and write to common transfer method in PatchGrid

  /**
    * Data access from all donor patches via m_InterCoeffData_WS
    * NOTE: Non padded version employing WeightedSet pattern.
    * ATTENTION: DOES NOT TURN ANY VECTORIAL VARIABLES !!
    * @param field the field, for which all variables are transfered
    */
  void accessDonorDataPadded(const size_t& field);


  /**
    * Data access from all donor patches via m_InterCoeffData_WS
    * Includes turning of one vectorial variable (example: speed vector) given by indicees (i_vx, i_vy, i_vz) in var type sequence
    * NOTE: donor_data and receive_data are small vectors with one entry per variable type to transfer
    *
    * @todo must improve variable selection for transformFree operation to allow better optimization.
    *
    * @param field the field, for which all variables are transfered
    * @param i_vx variable index forming the x-comp. of a vectorial variable. Example: "1" for speed (p,u,v,w,T)
    * @param i_vx variable index forming the y-comp. of a vectorial variable. Example: "2" for speed (p,u,v,w,T)
    * @param i_vx variable index forming the z-comp. of a vectorial variable. Example: "3" for speed (p,u,v,w,T)
    */
  void accessTurnDonorData_WS(const size_t& field,
                              const size_t& i_vx, const size_t& i_vy, const size_t& i_vz);


  void setFieldToZero(real *field);
  void setFieldToZero(size_t i_field);

  void setFieldToConst(real *field, real *var);
  void setFieldToConst(size_t i_field, real *var);

  void copyField(real *src, real *dst);
  void copyField(size_t i_src, size_t i_dst);
  void addField(real *src, real factor, real *dst);
  void addField(real *op1, real factor, real *op2, real *dst);
  void addField(size_t i_src, real factor, size_t i_dst);
  void addField(size_t i_op1, real factor, size_t i_op2, size_t i_dst);

  Transformation getTransformation() { return m_Transformation; }
  CoordTransformVV getTransformInertial2This() { return m_TransformInertial2This; }


  /** Compute ...
    * @return smallest characteristic length.
    */
  virtual real computeMinChLength() { BUG; return 0; }

  /**
   * @brief check if the GPU data pointer has been set already
   * @return true if the GPU data pointer has been set already
   */
  bool gpuDataSet() { return m_GpuDataSet; }

  /**
   * @brief get the GPU data pointer for pointer translation
   * @return the pointer to the data block on the GPU
   */
  real* getGpuData();

  /**
   * @brief ///< get indices of vectorial variables (e.g. veclovity, momentum, ...)
   * @return a vector with all starting indices
   */
  vector<size_t> getVectorVarIndices() { return m_VectorVarIndices; }

  /**
   * @brief Copy one field from host memory to GPU (device) memory
   * @param i_field the index of the field to copy
   */
  void copyFieldToDevice(size_t i_field);

  /**
   * @brief Copy one field from GPU (device) memory to host memory
   * @param i_field the index of the field to copy
   */
  void copyFieldToHost(size_t i_field);


  /**
   * Collect neighbour cells of a cell. Version for neighbour cells, that share a
   * common face with cell.
   * @param l_cell index of the cell for which to search neighbour cells.
   * @param l_cell_neighbours vector with neighbour cells indices (ret. ref.)
   */
  virtual void cellOverFaceNeighbours(const size_t& l_c,
                                      vector<size_t>& l_cell_neighbours) = 0;

  /**
   * Collect neighbour cells of a cell. Version for neighbour cells, that share at
   * least one common node with cell.
   * @param l_cell index of the cell for which to search neighbour cells.
   * @param l_cell_neighbours vector of neighbour cells of l_cell (ret. ref.)
   */
  virtual void cellOverNodeNeighbours(const size_t& l_c,
                                      vector<size_t>& l_cell_neighbours) = 0;


  /// @todo destructor
  virtual ~Patch();

};

inline void Patch::setFieldToZero(real *field)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    field[i] = 0.0;
  }
}

inline void Patch::setFieldToZero(size_t i_field)
{
  setFieldToZero(getField(i_field));
}

inline void Patch::setFieldToConst(real *field, real* var)
{
  for (size_t i_var = 0; i_var < m_NumVariables; ++i_var) {
    real* start_p = field + i_var*m_VariableSize;
    for (size_t i = 0; i < m_VariableSize; ++i) {
      start_p[i] = var[i_var];
    }
  }
  // data alignment error!!!
  // for (size_t i = 0; i < m_FieldSize; i=i+m_NumVariables) {
  //    for(size_t i_var = 0; i_var < m_NumVariables; ++i_var) {
  //      field[i+i_var] = var[i_var];
  //    }
  //  }
}


inline void Patch::setFieldToConst(size_t i_field, real *var)
{
  setFieldToConst(getField(i_field), var);
}


inline void Patch::copyField(real *src, real *dst)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    dst[i] = src[i];
  }
}


inline void Patch::copyField(size_t i_src, size_t i_dst)
{
  copyField(getField(i_src), getField(i_dst));
}


inline void Patch::addField(real *src, real factor, real *dst)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    dst[i] += factor*src[i];
  }
}


inline void Patch::addField(real *op1, real factor, real *op2, real *dst)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    dst[i] = op1[i] + factor*op2[i];
  }
}


inline void Patch::addField(size_t i_src, real factor, size_t i_dst)
{
  addField(getField(i_src), factor, getField(i_dst));
}


inline void Patch::addField(size_t i_op1, real factor, size_t i_op2, size_t i_dst)
{
  addField(getField(i_op1), factor, getField(i_op2), getField(i_dst));
}


inline vec3_t Patch::accessBBoxXYZoMin()
{
  if(!m_BBoxOk) {
    buildBoundingBox();
  }
  return m_BBoxXYZoMin;
}


/**
  * Access
  * @return upper coordinates of bounding box
  */
inline vec3_t Patch::accessBBoxXYZoMax()
{
  if(!m_BBoxOk) {
    buildBoundingBox();
  }
  return m_BBoxXYZoMax;
}

inline void Patch::copyFieldToDevice(size_t i_field)
{
#ifdef CUDA
  real *cpu_pointer = m_Data    + i_field*fieldSize();
  real *gpu_pointer = m_GpuData + i_field*fieldSize();
  cudaMemcpy(gpu_pointer, cpu_pointer, fieldSize()*sizeof(real), cudaMemcpyHostToDevice);
#endif
}

inline void Patch::copyFieldToHost(size_t i_field)
{
#ifdef CUDA
  real *cpu_pointer = m_Data    + i_field*fieldSize();
  real *gpu_pointer = m_GpuData + i_field*fieldSize();
  cudaMemcpy(cpu_pointer, gpu_pointer, fieldSize()*sizeof(real), cudaMemcpyDeviceToHost);
#endif
}


#endif // PATCH_H
