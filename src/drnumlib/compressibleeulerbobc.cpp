#include "compressibleeulerbobc.h"

CompressibleEulerBOBC::CompressibleEulerBOBC(size_t field)
  : BlockObjectBC(field)
{
}


CompressibleEulerBOBC::CompressibleEulerBOBC(size_t field, BlockObject* block_object)
  : BlockObjectBC(field, block_object)
{
}


void CompressibleEulerBOBC::operator()()
{
  BlockObject* bo = m_BlockObject;
  PatchGrid* pg = m_BlockObject->getPatchGrid();

  // Grey cells:
  // Indirect loop for patches affected (treats grey affection equally as black)
  for (size_t ii_p = 0; ii_p < bo->m_AffectedPatchIDs.size(); ii_p++) {
    Patch* patch = pg->getPatch(bo->m_AffectedPatchIDs[ii_p]);
    vector<greycell_t>& grey_cells = bo->m_GreyCells[ii_p];
    //.. Indirect loop for grey cells in patch
    for (size_t ll_cg = 0; ll_cg < grey_cells.size(); ll_cg++) {
      greycell_t& gc = grey_cells[ll_cg];
      size_t l_cell = gc.l_cell;
      //.... Set variable set to zero
      patch->setVarsetToZero (m_Field, l_cell);
      //.... loop for contributors (average of vars from sourrounding fluid nodes)
      for (size_t ll_cn = 0; ll_cn < gc.influencing_cells.size(); ll_cn++) {
        size_t l_cn = gc.influencing_cells[ll_cn];
        patch->addToVarset (m_Field, l_cell, l_cn);
      }
      patch->multVarsetScalar (m_Field, l_cell, gc.average_quote);
      //.... Set normal speed component to zero
      real varset_euler[5];
      patch->getVarset (m_Field, l_cell, varset_euler);
      real p, T, u, v, w;
      PerfectGas::conservativeToPrimitive(varset_euler, p, T, u, v, w);
      real n_speed = u * gc.n_vec[0] + v * gc.n_vec[1] + w * gc.n_vec[2];
      u -= n_speed * gc.n_vec[0];
      v -= n_speed * gc.n_vec[1];
      w -= n_speed * gc.n_vec[2];
      PerfectGas::primitiveToConservative(p, T, u, v, w, varset_euler);
      patch->setVarset (m_Field, l_cell, varset_euler);
    }
  }

  // Black cells:
  standardBlackCells ();

}
