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
#include "iteratorfeeder.h"


void IteratorFeeder::addIterator(PatchIterator* iterator)
{
  m_Iterators.push_back(iterator);
}

void IteratorFeeder::feed(PatchGrid &patch_grid)
{
  PatchGroups* patch_groups = patch_grid.getPatchGroups();

  // loop for patch groups
  for (size_t i_pg = 0 ; i_pg < patch_groups->accessNumPatchGroups(); ++i_pg) {
    SinglePatchGroup* spg = patch_groups->accessSinglePatchGroup(i_pg);
    CodeString cs = spg->m_SolverCodes;
    bool found = false;
    for (size_t i_it = 0; i_it < m_Iterators.size(); ++i_it) {
      /** @todo Direct comparison may void tolerance concerning keyword "0" .
        * This happens if not any of the patches has requested a non "0"
        * flux keyword at the corresponding position, while the iterator has been
        * build with it.
        * Example: iterator has CodeString = "a b c d"
        *          all patches request CodeStrig = "a b c 0"
        *  => no iterator will be found, causing "BUG;" below.
        * How to treat this?
        * Case could as well be considered an error, as an iterator has been
        * build with a flux that is never used.
        */
      if (m_Iterators[i_it]->getCodeString() == cs) {
        found = true;
        for (size_t i_patch = 0; i_patch < spg->m_Patches.size(); ++i_patch) {
          m_Iterators[i_it]->addPatch(spg->m_Patches[i_patch]);
        }
      }
    }
    if (!found) {
      cout << cs.c_str() << endl;
      BUG;
    }
  }
}
