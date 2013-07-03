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
      BUG;
    }
  }
}
