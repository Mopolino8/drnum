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
#include "compressiblesolidwallbobc.h"

CompressibleSolidWallBOBC::CompressibleSolidWallBOBC (size_t field)
  : BlockObjectBC (field)
{
}


CompressibleSolidWallBOBC::CompressibleSolidWallBOBC (size_t field, BlockObject* block_object)
  : BlockObjectBC(field, block_object)
{
}


void CompressibleSolidWallBOBC::operator()()
{
  BlockObject* bo = m_BlockObject;
  PatchGrid* pg = m_BlockObject->getPatchGrid();

  // Grey cells:
  // Indirect loop for patches affected (treats grey affection equally as black)
  for (size_t ii_p = 0; ii_p < bo->getAffectedPatchIDs().size(); ii_p++) {
    Patch* patch = pg->getPatch(bo->getAffectedPatchIDs()[ii_p]);
    vector<greycell_t>& grey_cells = bo->getGreyCells()[ii_p];
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
      //.... Extrapolate speed to get virtual "0" at presumed interface.
      //     See document xxxx?
      real varset[5];
      patch->getVarset (m_Field, l_cell, varset);
      real p, T, u, v, w;
      PerfectGas::conservativeToPrimitive(varset, p, T, u, v, w);
      //...... geometric quotient of extrapolation
      real gf = gc.grey_factor;
      real quot = (1. - 2.*gf) / (3. - 2.*gf);
      u *= quot;
      v *= quot;
      w *= quot;
      PerfectGas::primitiveToConservative(p, T, u, v, w, varset);
      patch->setVarset (m_Field, l_cell, varset);
    }
  }

  // Black cells:
  /// @todo A first order extrapolation of 1st layer black cells would be better.
  standardBlackCells0 ();

}
