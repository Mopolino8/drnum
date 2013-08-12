#include "blockobjectbc.h"

BlockObjectBC::BlockObjectBC (size_t field)
{
  m_Field = field;
}


BlockObjectBC::BlockObjectBC (size_t field, BlockObject* block_object)
{
  m_Field = field;
  setBlockObject (block_object);
}


void BlockObjectBC::standardBlackCells0 ()
{
  BlockObject* bo = m_BlockObject;
  PatchGrid* pg = m_BlockObject->getPatchGrid();

  // Loop for black groups
  for (size_t i_bg = 0; i_bg < bo->m_NumBlackGroups; i_bg++) {
    //.. Loop for patches
    for (size_t ii_p = 0; ii_p < bo->m_AffectedPatchIDs.size(); ii_p++) {
      Patch* patch = pg->getPatch(bo->m_AffectedPatchIDs[ii_p]);
      vector<blackcell_t>& black_cells = bo->m_BlackCells[i_bg][ii_p];
      //.. Indirect loop for black cells in patch
      for (size_t ll_cb = 0; ll_cb < black_cells.size(); ll_cb++) {
        blackcell_t& bc = black_cells[ll_cb];
        size_t l_cell = bc.l_cell;
        //.... Set variable set to zero
        patch->setVarsetToZero (m_Field, l_cell);
        //.... loop for contributors (average of vars from sourrounding fluid nodes)
        for (size_t ll_cn = 0; ll_cn < bc.influencing_cells.size(); ll_cn++) {
          size_t l_cn = bc.influencing_cells[ll_cn];
          patch->addToVarset (m_Field, l_cell, l_cn);
        }
        patch->multVarsetScalar (m_Field, l_cell, bc.average_quote);
      }
    }
  }
}

