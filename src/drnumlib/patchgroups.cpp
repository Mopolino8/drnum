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
//#include <unistd.h>
#include "patchgroups.h"


PatchGroups::PatchGroups(bool allow_solver_exceptions)
{
  m_AllowSolverExceptions = allow_solver_exceptions;
}


pair<size_t, size_t> PatchGroups::insertPatch(Patch* patch)
{
  return insertPatch(patch, m_AllowSolverExceptions);
}


pair<size_t, size_t> PatchGroups::insertPatch(Patch* patch, bool allow_solver_exceptions)
{
  size_t groupindex = 0;
  size_t patchindex = 0;
  size_t patchtype = patch->accessPatchType();
  // Check all previously created SinglePatchGroups aiming to find same patchtype and suitable
  // solvercodes.
  // Note performance: assume these are few.
  bool found = false;
  for(size_t i = 0; i < m_PatchGroups.size(); i++) {
    // Check, if patch type matches
    if(m_PatchGroups[i]->m_PatchType == patchtype) {
      // Check, if solvercodes match:
      //   allsame == true    : solver codes match exactly.
      //   overruled == true  : one or more solver codes of patch::m_SolverCodes are "0",
      //                        while corresponding ones in m_PatchGroups[i].m_SolverCodes arent.
      //   underruled == true : one or more solver codes in m_PatchGroups[i].m_SolverCodes are
      //                        "0", while the corresponding in patch::m_SolverCodes arent.
      //   concatsolvercodes  : concatenated solver codes, replacing the "zeros" by
      //                        living codes
      bool allsame, overruled, underruled;
      CodeString concatsolvercodes;
      CodeString patchsolvercodes = patch->accessSolverCodes();
      concatsolvercodes = patchsolvercodes.compareCodeString(m_PatchGroups[i]->m_SolverCodes,
                                                             allsame,
                                                             overruled,
                                                             underruled);

//      patch->compareSolverCodes(m_PatchGroups[i]->m_SolverCodes,
//                                allsame,
//                                overruled,
//                                underruled;
//                                concatsolvercodes);
      // Check, if SinglePatchGroup can be used
      bool usegroup = allsame;
      usegroup = usegroup || ((overruled || underruled) && allow_solver_exceptions);
      // Insert, if useful
      if(usegroup) {
        found = true;
        groupindex = i;
        patchindex = m_PatchGroups[i]->insert(patch);
        //.. Adapt solver codes, if needed
        if(underruled) {
          m_PatchGroups[i]->m_SolverCodes = concatsolvercodes;
        }
        if(underruled || overruled) {
          m_PatchGroups[i]->m_SolverExceptions = true;
        }
        //.. found suitable group, stop search
        break;
      }
    }
  }
  // Create a new SinglePatchGroup, if none of the previous suits
  if(!found) {
    // create a new SinglePatchGroup
    SinglePatchGroup* newspg;
    newspg = new SinglePatchGroup();
    newspg->m_PatchType = patchtype;
    newspg->m_SolverCodes = patch->accessSolverCodes();
    newspg->m_SolverExceptions = false;
    patchindex = newspg->insert(patch);
    m_PatchGroups.push_back(newspg);
    groupindex = m_PatchGroups.size() - 1;
  }
  // return value
  pair<size_t, size_t> retval;
  retval.first = groupindex;
  retval.second = patchindex;
  return retval;
}

void PatchGroups::insertPatchGrid(PatchGrid* patchgrid, bool allow_solver_exceptions)
{
  size_t numpatches = patchgrid->getNumPatches();
  for(size_t ip = 0; ip < numpatches; ip++) {
    insertPatch(patchgrid->getPatch(ip), allow_solver_exceptions);
  }
}

void PatchGroups::insertPatchGrids(vector<PatchGrid*> patchgrids, bool allow_solver_exceptions)
{
  for(size_t ipg = 0; ipg < patchgrids.size(); ipg++) {
    insertPatchGrid(patchgrids[ipg], allow_solver_exceptions);
  }
}

PatchGroups::~PatchGroups()
{
  /// @todo missing
}
