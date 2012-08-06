#ifndef PATCHGRID_H
#define PATCHGRID_H

#include <cstddef>

#include "blockcfd.h"
#include "TList.hh"
#include "TInsectionList.hh"
#include "patch.h"
/// @todo thought I have it from oversetmouse, but didnt find it yet, other name? => write new
//#include "hashraster.h"

class PatchGrid
{

protected: // attributes

  /** List of patches in the grid */
  TList<Patch*>* m_patches;

  /// @todo might eliminate this info and store in Patch class
  /** Grid dependencies. These are the other patches that contribute
   *  to the interpolation sets. */
  TInsectionList<Patch*>* m_patchdependencies;

  vec3_t m_bbox_xyzo_min;   ///< lowest coordinates of smallest box around grid in inertial coords.
  vec3_t m_bbox_xyzo_max;   ///< highest coordinates of smallest box around grid in inertial coords.
  bool m_bbox_OK;           ///< flag indicating wether the bounding box is available
  bool m_dependencies_OK;   ///< indicate, if dependencies are updated

private: // methods

protected: // methods

public: // methods

  PatchGrid();

  /**
   * Initialize patch lists
   * @param size_default default patch list size (number of patches to hold)
   * @param size_incr list rowth increment
   */
  void initLists(size_t size_default, size_t size_incr);

  /** Insert a new patch
   * @param new_patch patch to insert
   * @return the index of the patch in the m_patchs
   */
  size_t insertPatch(Patch* new_patch)
  {
    size_t i_patch = m_patches->AddEntry();
    m_patches->At(i_patch) = new_patch;
    /// @todo how about dependencies??
    m_dependencies_OK = false;
    return i_patch;
  }

  /** Delete a patch
   * @param i_patch the index of the patch in m_patches
   */
  void deletePatch(size_t i_patch)
  {
    m_patches->DelEntry(i_patch);
    m_dependencies_OK = false;
  }

  /**
   * Read patch list from file.
   */
  void readGrid() {BUG;}

  /**
   * Write patch list to file.
   */
  void writeGrid() {BUG;}

  /**
   * Find the patch dependencies, if not given as input.
   * @param with_intercoeff flag directive to build transfer coefficients
   */
  void computeDependencies(bool with_intercoeff);

  /**
    * Build bounding box around whole grid
    */
  void buildBoundingBox();

  // Access methods
  size_t getNumPatches() {return m_patches->NumEntries();}
  TList<Patch*>* getPatches() {return m_patches;}
  TInsectionList<Patch*>*getPatchdependencies() {return m_patchdependencies;}

  virtual ~PatchGrid();

};

#endif // PATCHGRID_H
