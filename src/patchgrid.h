#ifndef PATCHGRID_H
#define PATCHGRID_H

#include <cstddef>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

class PatchGrid;

#include "blockcfd.h"
#include "patch.h"
#include "cartesianpatch.h"
#include "vectorhashraster.h"
#include "patchgroups.h"

class PatchGrid
{

protected: // attributes

  vector<Patch*> m_patches;                ///< List of patches in the grid
  PatchGroups* m_patchgroups;

  //VectorHashRaster<size_t> m_HashRaster;   ///< Hash raster to assist orientation

  // settings (same as individually defined in patch.h)

 /// @todo m_NumFields m_NumVariables must go ...
  size_t  m_NumFields;       ///< number of fields (e.g. old, new, ...)
  size_t  m_NumVariables;    ///< number of variables (e.g. rho, rhou, ...)
  bool m_InterpolateData;    ///< Flag indicates wether to interpolate data on interpatch transfers
  bool m_InterpolateGrad1N;  ///< Flag indicates wether to interpolate directed gradients on interpatch transfers
  bool m_TransferPadded;     ///< Flag indicates wether to transfer donor data in padded versions with "InterCoeffPad".
  size_t m_NumProtectLayers; ///< number of boundary protection layers, in which no interpol access from other patches is allowed
  size_t m_NumOverlapLayers; ///< number of boundary cell layers, for which to get data from donor neighbour patches

  vec3_t m_bbox_xyzo_min;   ///< lowest coordinates of smallest box around grid in inertial coords.
  vec3_t m_bbox_xyzo_max;   ///< highest coordinates of smallest box around grid in inertial coords.
  bool m_bbox_OK;           ///< flag indicating wether the bounding box is available
  bool m_dependencies_OK;   ///< indicate, if dependencies are updated

private: // methods

protected: // methods

public: // methods

  PatchGrid(size_t num_protectlayers=1, size_t num_overlaplayers=1);

  /** Set number of variable fields on all patches of patchgrid.
    * @param num_fields number of variable fields.
    */
  void  setNumberOfFields(size_t num_fields);

  /**
    * Set number of variables on all patches of patchgrid.
    * @param num_variables number of variables.
    */
  void  setNumberOfVariables(size_t num_variables);

  /**
    * Set interaction with/without data transfers
    * @param interpolate_data bool to cause data interpolation if true
    */
  void setInterpolateData(bool interpolatedata = true);

  /**
    * Set interaction with/without 1. gradient transfers
    * @param interpolate_data bool to cause gradient interpolation if true
    */
  void setInterpolateGrad1N(bool interpolategrad1N = true);

  /**
    * Set all dependency transfers from any donors to be padded, employing data transfer
    * classes of type "InterCoeffPad".
    * @param trans_padded bool to cause padded data transfers
    */
  void setTransferPadded(bool trans_padded = true);

  /**
    * Set number of protection layers
    * @param num_protectlayers number of protection layers
    */
  void setNumProtectLayers(size_t num_protectlayers);

  /**
    * Set number of protection layers
    * @param num_overlaplayers number of overlap layers
    */
  void setNumOverlapLayers(size_t num_overlaplayers);

  /** Insert a new patch
   * @param new_patch patch to insert
   * @return the index of the patch in the m_patchs
   */
  size_t insertPatch(Patch* new_patch);

  /**
    * Hand over general attributes (logical settings) to a patch
    */
  void setGeneralAttributes(Patch* patch);

  /** Delete a patch
   * @param i_patch the index of the patch in m_patches
   */
  void deletePatch(size_t i_patch);

  /**
   * Read patch list from file.
   * @param gridfilename filename of grid file relative to cwd.
   */
  void readGrid(string gridfilename = "/grid/patches");

  /**
   * Write patch list to file.
   */
  void writeGrid() {BUG;}

  /**
    * Scale total patchgrid relative to origin of parental coordsyst.
    * NOTE: Affects reference positions and physical patch sizes.
    * @param scfactor scaling factor.
    */
  void scaleRefParental(real scfactor);

  /**
   * Find the patch dependencies, if not given as input.
   * @param with_intercoeff flag directive to build transfer coefficients
   */
  void computeDependencies(const bool& with_intercoeff);

  /**
     * Envoque same function for all patches.
     * - Reduce contribution weights for receiving cells, being influenced by more than one donor patch.
     * - Transfer data to padded data sets, if required.
     */
  void finalizeDependencies();

  /**
    * Envoque data access to neighbour patches for all patches in the grid.
    * Version using m_InterCoeffData_WS data
    * NOTE: Non padded version employing WeightedSet pattern.
    * ATTENTION: DOES NOT TURN ANY VECTORIAL VARIABLES !!
    * @param field the field, for which all variables are transfered
    */
  void accessAllDonorData_WS(const size_t& field);


  /** @todo Let patch grid know, if any vectorial variables must be turned (default:yes) and which variables are affected
    *       (idicees in variable set, like 1,2,3 for speeds in (p, u, v, w, ...) */


  /**
    * Build bounding box around whole grid
    */
  void buildBoundingBox(const bool& force = true);

  /**
    * Build a hash box around whole grid
    */
  //void buildHashRaster(size_t resolution = 1000000, bool force = true);
  void buildHashRaster(size_t resolution, bool force,
                       VectorHashRaster<size_t>& m_HashRaster);
  /**
    * Initialize. variables: Copy a variable set onto each cell of all patches.
    * @param i_field variable field on which to copy data
    * @param var single set of variables, example (u, v, w, p).
    */
  void setFieldToConst(size_t i_field, real *var);

  // Access methods

  size_t getNumPatches() {return m_patches.size();}

  Patch* getPatch(const size_t& ip) {return m_patches[ip];}

  PatchGroups* getPatchGroups() {return m_patchgroups;}

  SinglePatchGroup* getSinglePatchGroup(const size_t& ipg);

  /**
    * Compute ...
    * @return smallest length in any of the patches.
    */
  real computeMinChLength();

  virtual ~PatchGrid();

};

#endif // PATCHGRID_H
