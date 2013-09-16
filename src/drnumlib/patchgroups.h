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
#ifndef PATCHGROUPS_H
#define PATCHGROUPS_H

#include <cstddef>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

struct SinglePatchGroup;
class PatchGroups;

#include "blockcfd.h"
#include "patch.h"
#include "patchgrid.h"
#include "codestring.h"

/**
  * SinglePatchGroup characterises a group of patches with following (mostly) common properties:
  *   1) Geometric Type (m_PatchType): example: all of type CartesianPatch
  *   2) Solver code string:
  *      All patches obei same computations as given by code words in m_SolverCodes
  *      Exceptions for patches might be allowed, replacing one or more code words
  *      by "0". If code words are "0", this is equivalent to flagging out the corresponding
  *      action.
  */

/// @todo make class with access methods
struct SinglePatchGroup
{
  size_t m_PatchType;       ///< Geometric type code of all patches in patch-group
  vector<Patch*> m_Patches; ///< List of patches-indicees in sequence of PatchGrid::m_Patches
  CodeString m_SolverCodes; ///< Code-string for the patch group
  bool m_SolverExceptions;  ///< True if any solver exceptions occured on this group yet

  /**
    * Insert
    * @param patch ptr of patch
    * @return index of patch in m_Patches
    */
  size_t insert(Patch* patch)
  {
    m_Patches.push_back(patch);
    return m_Patches.size() - 1;
  }

};


/**
  *  PatchGroups holds all SinglePatchGroup-structs.
  */
class PatchGroups
{

protected: // attributes

  vector<SinglePatchGroup*> m_PatchGroups; ///< List of patch-groups
  bool m_AllowSolverExceptions;    ///< True: solver exceptions are generally allowed unless overwritten

private: // methods

protected: // methods

public: // methods

  /**
    * Empty constructor.
    */
  PatchGroups(bool allow_solver_exceptions = true);

  /** Insert a patch. Will look for a matching SinglePatchGroup, previously available and
    * insert there. If none is found a new SinglePatchGroup is created.
    * Note indexing: groupindex: index to m_PatchGroups
    *                patchindex in group: index to SinglePatchGroup::m_Patches
    * @param patch patch to insert
    * @param allowsolverexceptions indicates that code exceptions are allowed, see SinglePatchGroup.
    * @return pair of indixees (groupindex, patchindex in group)
    */
  pair<size_t, size_t> insertPatch(Patch* patch, bool allow_solver_exceptions);

  /** Insert a patch. Will look for a matching SinglePatchGroup, previously available and
    * insert there. If none is found a new SinglePatchGroup is created.
    * @param patch patch to insert
    * @param groupindex index of SinglePatchGroup in list m_PatchGroups (return reference)
    * @param patchindex within group, sequence of SinglePatchGroup::m_Patches (return reference)
    */
  pair<size_t, size_t> insertPatch(Patch* patch);

  /** Build up from an entire PatchGrid
    * @param patchgrid the PatchGrid from which to build
    * @param allowsolverexceptions indicates that bc exceptions are allowed
    */
  void insertPatchGrid(PatchGrid* patchgrid, bool allow_solver_exceptions);

  /** Build up from a vector of PatchGrid. Needed for overset computations.
    * @param patchgrids the vector of PatchGrid from which to build
    * @param allowsolverexceptions indicates that bc exceptions are allowed
    */
  void insertPatchGrids(vector<PatchGrid*> patchgrids, bool allow_solver_exceptions);

  /** Access
    * @return number of patch groups
    */
  size_t accessNumPatchGroups()
  {
    return m_PatchGroups.size();
  }

  /** Access a SinglePatchGroup
    * @param spg_index index of SinglePatchGroup
    * @return ptr to SinglePatchGroup
    */
  SinglePatchGroup* accessSinglePatchGroup(const size_t& spg_index)
  {
    return m_PatchGroups[spg_index];
  }

  /** Accessnumber of patch groups of a given geometric type
    * @param patchtype patch type
    * @return number of patch groups of a given geometric type
    */
  size_t accessNumPatchGroups(size_t patch_type)
  {
    size_t numspg = 0;
    for(size_t i = 0; i < m_PatchGroups.size(); i++)
    {
      if(m_PatchGroups[i]->m_PatchType == patch_type) {
        numspg++;
      }
    }
    return numspg;
  }

  /** Write to file
    * @param filename filename to write to
    */
  void writeToFile(string filename);


  virtual ~PatchGroups();

};

#endif // PATCHGROUPS_H
