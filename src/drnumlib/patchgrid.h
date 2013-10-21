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
#ifndef PATCHGRID_H
#define PATCHGRID_H

#include <cstddef>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

class PatchGrid;

#include "drnum.h"
#include "patch.h"
#include "cartesianpatch.h"
#include "vectorhashraster.h"
#include "patchgroups.h"
#include "postprocessingvariables.h"


/// @todo uncommented Grad1N stuff. Clean up.


class PatchGrid
{

protected: // attributes

  vector<Patch*> m_Patches;      ///< List of patches in the grid
  PatchGroups* m_PatchGroups;

 /// @todo m_NumFields m_NumVariables must go ...
  size_t m_NumFields;           ///< number of fields (e.g. old, new, ...)
  size_t m_NumVariables;        ///< number of variables (e.g. rho, rhou, ...)
  bool   m_InterpolateData;     ///< Flag indicates wether to interpolate data on interpatch transfers
  bool   m_TransferPadded;      ///< Flag indicates wether to transfer donor data in padded versions with "InterCoeffPad".
  string m_TransferType;        ///< Keyword indicating transfer type ("ws", "padded", "padded_direct")
  size_t m_NumSeekLayers;       ///< number of boundary cell layers, for which to get data from donor neighbour patches
  size_t m_NumAddProtectLayers; ///< additional number of boundary protected layers, in which no interpol access from other patches is allowed

  vec3_t m_BboxXyzoMin;         ///< lowest coordinates of smallest box around grid in inertial coords.
  vec3_t m_BboxXyzoMax;         ///< highest coordinates of smallest box around grid in inertial coords.
  bool   m_BboxOk;              ///< flag indicating wether the bounding box is available
  bool   m_DependenciesOk;      ///< indicate, if dependencies are updated

  vector<size_t> m_VectorVarIndices;  ///< definition of vectorial variables by index corresponding x-direction. Subsequent two: y and z

private: // methods

protected: // methods

public: // methods

  /** Constructor
    * @param num_seeklayers default number of element layers on patch boundaries,
    *        seekinng interpolated data
    * @param num_addprotectlayers default number of additional, protected element layers
    *        between seek layer and core region
    */
  PatchGrid(size_t num_seeklayers = 2, size_t num_addprotectlayers = 0);


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
    * Define a variable triple as vectorial variable.
    * This is done by setting the first index (direction: local x) of the vectorial variable.
    * Subsequent two variables correspond to local y and z direcions.
    * @param index_x the index of x-coord of vectorial variable
    */
  void defineVectorVar(const size_t& index_x);


  /**
    * Set all dependency transfers to a type.
    * Known types are "ws" (weightedset), "padded" (InterCoeffPad) and
    * "padded_direct" (use direct padded lists
    * @param trans_type string indicating transfer mode
    */
  void setTransferType(string trans_type = "padded_direct");


  /// @todo The following 2 methods might perhaps go, as default is set construcion time (?)
  /**
    * Set default number of additional protection layers
    * @param num_protectlayers number of protection layers
    */
  void setNumAddProtectLayers(size_t num_addprotectlayers);


  /**
    * Set default number of element layers, seeking data from other patches.
    * @param num_seeklayers number of seeking element layers
    */
  void setNumSeekLayers(size_t num_seeklayers);


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
  //void deletePatch(size_t i_patch);
  void deletePatch(size_t) {
    /** @todo Other patches may be dependent. These dependencies must be cleared.
      * Most likely this method will never be required, unless for mesh generation purposes.
      */
    BUG;
  }


  /**
   * Read patch list from file.
   * @param gridfilename filename of grid file relative to cwd.
   */
  void readGrid(string gridfilename = "/grid/patches", real scale = 1.0);


  /// @todo not implemented
  /**
   * Write patch list to file.
   */
  void writeGrid() {BUG;}


  /**
   * This has been changed in order to have all patches in a single file.
   * The last commit before the change is: c62f01bee34b3832b648a133353241740ea7a835
   * @param i_field the field index of the field to be read
   * @param file_name full file name relative to cwd.
   * @return the last simulation time
   */
  real readData(size_t i_field, QString file_name);


  /**
   * This has been changed in order to have all patches in a single file.
   * The last commit before the change is: c62f01bee34b3832b648a133353241740ea7a835
   * @param i_field the field index of the field to be written
   * @param base_file_name file name relative to cwd.
   * @param time the current simulation time
   * @param discrete counter (usually time counter).
   */
  void writeData(size_t i_field, QString base_file_name, real time, int count);


  /**
   * @brief Write all patches to a VTK multi-block file.
   * @param file_name the file name ("*.vtm" will be added)
   * @param proc_vars this defines which varaibles will be written for post-processing
   */
  void writeToVtk(size_t i_field, string file_name, const PostProcessingVariables &proc_vars, int count = -1);


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
    * Write a log file for patch dependencies.
    * @param dep_log_filename dependencies log filename
    */
  void writeDependenciesLog(string dep_log_filename = "./patches/dependencies_log");


  /**
    * Find list of patches, that overlap a given cartesian box in
    * xyzo-system.
    * @param cbox_min lower coords of box
    * @param cbox_min upper coords of box
    * @param only_core indicates to analise only core region of patch
    * @param overlap_patches vector of patch indices overlapping box (ret. ref.)
    */
  void findBoxOverlappingPatches(const vec3_t& cbox_min, const vec3_t& cbox_max,
                                 const bool& only_core,
                                 vector<size_t>& overlap_patches);

  /// @todo Change to virtual function or type sorted data transfers.

  /**
    * Envoque data access to neighbour patches for all patches in the grid.
    * Does own transfer type selection upon m_TransferType .
    * ATTENTION: DOES NOT TURN ANY VECTORIAL VARIABLES !!
    * @param field the field, for which all variables are transfered
    */
  void accessAllDonorData(const size_t& field);


  /**
    * Envoque data access to neighbour patches for all patches in the grid.
    * Version using m_InterCoeffData_WS data
    * NOTE: Non padded version employing WeightedSet pattern.
    * ATTENTION: DOES NOT TURN ANY VECTORIAL VARIABLES !!
    * @param field the field, for which all variables are transfered
    */
  void accessAllDonorData_WS(const size_t& field);


  /**
    * Envoque data access to neighbour patches for all patches in the grid.
    * Version using m_InterCoeffData data
    * NOTE: Non padded version employing WeightedSet pattern.
    * ATTENTION: DOES NOT TURN ANY VECTORIAL VARIABLES !!
    * @param field the field, for which all variables are transfered
    * @param direct bool, indicating wether to use direct transfer lists
    */
  void accessAllDonorDataPadded(const size_t& field, const bool& direct);


  /** @todo Let patch grid know, if any vectorial variables must be turned (default:yes) and which variables are affected
    *       (idicees in variable set, like 1,2,3 for speeds in (p, u, v, w, ...) */


  /**
    * Build bounding box around whole grid
    */
  void buildBoundingBox(const bool& force = true);


  /**
    * Build a hash box around whole grid
    */
  void buildHashRaster(size_t resolution, bool force,
                       VectorHashRaster<size_t>& m_HashRaster);


  /**
    * Initialize. variables: Copy a variable set onto each cell of all patches.
    * @param i_field variable field on which to copy data
    * @param var single set of variables, example (u, v, w, p).
    */
  void setFieldToConst(size_t i_field, real *var);


  // Access methods

  size_t getNumPatches() {return m_Patches.size();}


  Patch* getPatch(const size_t& ip) {return m_Patches[ip];}


  PatchGroups* getPatchGroups() {return m_PatchGroups;}


  SinglePatchGroup* getSinglePatchGroup(const size_t& ipg);


  /**
    * Compute ...
    * @return smallest length in any of the patches.
    */
  real computeMinChLength();

  vector<size_t> getVectorVarIndices();

  bool findCell(vec3_t xo, int &id_patch, int &id_cell);

  virtual ~PatchGrid();

};

#endif // PATCHGRID_H
