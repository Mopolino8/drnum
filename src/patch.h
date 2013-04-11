#ifndef PATCH_H
#define PATCH_H

//#include <cstddef>
//#include <string.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>

#include "blockcfd.h"

#include <iostream>
#include <sstream>
#include <string>
using namespace std;

//#include "blockcfd.h"

#include "utility/weightedset.h"

class Patch;

#include "intercoeffpad.h"
#include "intercoeffws.h"
#include "math/coordtransformvv.h"

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

#include "patch_common.h"

  void  allocateData();

protected: // attributes

  size_t m_mytypecode;        ///< Type-code of a patch, e.g. 1001: CartesianPatch, etc...
  string m_solvercodes;       ///< String containing infos for choice of field-fluxes, RB-fluxes, sources, ...

  // reference position
  /// @todo redundant info: ommit m_Xo,m_Yo,m_Zo as it is contained in m_transformInertial2This
  real m_Xo;                  ///< X-coord of reference position in parental coordinates
  real m_Yo;                  ///< Y-coord of reference position in parental coordinates
  real m_Zo;                  ///< Z-coord of reference position in parental coordinates
  // settings
  bool m_InterpolateData;     ///< Flag indicates wether to interpolate data on interpatch transfers
  bool m_InterpolateGrad1N;   ///< Flag indicates wether to interpolate directed gradients on interpatch transfers
  bool m_TransferPadded;      ///< Flag indicates wether to transfer donor data in padded versions with "InterCoeffPad".
  bool m_ProtectException ;   ///< Flag indicates a protection exception: the number of protection layers is not constant.
  size_t m_NumProtectLayers;  ///< number of boundary protection layers, in which no interpol access from other patches is allowed
  size_t m_NumOverlapLayers;  ///< number of boundary cell layers, for which to get data from donor neighbour patches

  // IO-scaling, position, orientation
  real m_ioscale;                            ///< general scale factor of mesh related to io-values
  CoordTransformVV m_transformInertial2This; ///< transformation matrix to transform intertial coords into system of "this"

  // bounding box
  vec3_t m_bbox_xyzo_min;                    ///< lowest coordinates of smallest box around patch in inertial coords.
  vec3_t m_bbox_xyzo_max;                    ///< highest coordinates of smallest box around patch in inertial coords.
  bool m_bbox_OK;                            ///< flag indicating wether the bounding box is available

  // intermediate variables
  bool m_receiveCells_OK; ///< Flag indicating that receive_cells have been extracted yet. Might be set false on mesh changes.

  // lists related to receiving cells in overlap layers
  vector<size_t> m_receive_cells;           ///< cells of "this", expecting to get data (and/or grads) from any donor neighbour
  vector<size_t> m_receive_cell_data_hits;  ///< number of contributing patches for data, note indexing as receive_cells
  vector<size_t> m_receive_cell_grad1N_hits;///< number of contributing patches for grad1N, note indexing as receive_cells

  // lists related to neighbouring donor patches
  vector<pair<Patch*, CoordTransformVV> > m_neighbours; ///< neighbouring donor patches and coord transformation. NOTE: receiving patches not stored.
  vector<InterCoeffPad> m_InterCoeffData;               ///< Interpolation coefficient lists for data
  vector<InterCoeffPad> m_InterCoeffGrad1N;             ///< Interpolation coefficient lists for 1st directed gradients
  vector<InterCoeffWS> m_InterCoeffData_WS;             ///< same as m_InterCoeffData, but with WeightedSets (CPU only)
  vector<InterCoeffWS> m_InterCoeffGrad1N_WS;           ///< same as m_InterCoeffGrad1N, but with WeightedSets (CPU only)
  // postponed  vector<InterCoeff*> m_InterCoeffGrad2N;    ///< Interpolation coefficient lists for 2nd directed gradients


protected: // methods

  // internal data handling
  void  deleteData();
  void  resize(size_t variable_size);

  void setTransformation(Transformation t) { m_Transformation = t; }  /// @todo keep for compatibility, prefer CoordTransformVV later

  // geometry
  virtual void buildBoundingBox()=0;

  // external data handling
  virtual void extractReceiveCells() { BUG; } ///< Extract cells in overlap layers (m_receive_cells)
  void compactReceiveCellLists();             ///< Clean up m_receive_cells and hit counter lists


public: // methods

  /**
   * Constructor
   * NOTE: The default number of protection as well as overplap layers is 1
   * @param num_protectlayers number of protection layers, in which no foreign data access is allowed
   * @param num_overlaplayers number of overlap layers, for which to get data from donor neighbour patches
   */
  Patch(size_t num_protectlayers=1, size_t num_overlaplayers=1);

  virtual ~Patch();

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

  /**
    * Set number of protection layers
    * @param num_protectlayers number of protection layers
    */
  virtual void setNumProtectLayers(size_t num_protectlayers)
  {
    m_NumProtectLayers = num_protectlayers;
  }

  /**
    * Set number of protection layers
    * @param num_overlaplayers number of overlap layers
    */
  void setNumOverlapLayers(size_t num_overlaplayers)
  {
    m_NumOverlapLayers = num_overlaplayers;
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
    * Set interaction with/without 1. gradient transfers
    * @param interpolate_data bool to cause gradient interpolation if true
    */
  void setInterpolateGrad1N(bool interpolategrad1N = true)
  {
    m_InterpolateGrad1N = interpolategrad1N;
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
    * @param s_mesh the stream to read from
    * @return true, if successful
    */
  virtual bool readFromFile(istringstream iss_input);

  /**
    * Write mesh data to file
    * @param s_mesh the stream to write to
    * @return true, if successful
    */
  virtual bool writeToFile(ifstream&) {return true;};
//  virtual bool writeToFile(ifstream& s_mesh) {return true;};

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
    m_transformInertial2This.scaleVector(scfactor);
  }

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
  for (size_t i = 0; i < m_FieldSize; i=i+m_NumVariables) {
    for(size_t i_var = 0; i_var < m_NumVariables; ++i_var) {
      field[i+i_var] = var[i_var];
    }
  }
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
  if(!m_bbox_OK) {
    buildBoundingBox();
  }
  return m_bbox_xyzo_min;
}

/**
  * Access
  * @return upper coordinates of bounding box
  */
inline vec3_t Patch::accessBBoxXYZoMax()
{
  if(!m_bbox_OK) {
    buildBoundingBox();
  }
  return m_bbox_xyzo_max;
}

#endif // PATCH_H
