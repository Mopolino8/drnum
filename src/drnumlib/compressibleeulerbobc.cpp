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
#include "compressibleeulerbobc.h"

CompressibleEulerBOBC::CompressibleEulerBOBC (size_t field)
  : BlockObjectBC (field)
{
}


CompressibleEulerBOBC::CompressibleEulerBOBC (size_t field, BlockObject* block_object)
  : BlockObjectBC (field, block_object)
{
}


void CompressibleEulerBOBC::operator()()
{
  // ATTENTION Assume first 5 variables to make up compressible vars.
  /// @todo Will overwrite variables beyond the 5th with undefined code

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
      patch->setVarsubsetToZero (m_Field, l_cell,
                                 0, 5);
      //.... loop for contributors (average of vars from sourrounding fluid nodes)
      for (size_t ll_cn = 0; ll_cn < gc.influencing_cells.size(); ll_cn++) {
        size_t l_cn = gc.influencing_cells[ll_cn];
        patch->addToVarsubset (m_Field, l_cell, l_cn,
                               0, 5);
      }
      patch->multVarsubsetScalar (m_Field, l_cell, gc.average_quote,
                                  0, 5);
      //.... Set normal speed component to zero
      real varset_euler[5];
      patch->getVarsubset (m_Field, l_cell, varset_euler,
                           0, 5);
      real p, T, u, v, w;
      PerfectGas::conservativeToPrimitive(varset_euler, p, T, u, v, w);
      real n_speed = u * gc.n_vec[0] + v * gc.n_vec[1] + w * gc.n_vec[2];
      u -= n_speed * gc.n_vec[0];
      v -= n_speed * gc.n_vec[1];
      w -= n_speed * gc.n_vec[2];
      PerfectGas::primitiveToConservative(p, T, u, v, w, varset_euler);
      patch->setVarsubset (m_Field, l_cell, varset_euler,
                           0, 5);
    }
  }

  // Black cells:
  standardBlackCells0 ();

}
