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
#ifndef CPU_LEVELSETOBJECTBC_H
#define CPU_LEVELSETOBJECTBC_H

#include "levelsetobjectbc.h"
#include "lslayerdataextrapol.h"

template <typename OP>
class CPU_LevelSetObjectBC : public LevelSetObjectBC
{

protected:
  OP m_Op;

public:

  CPU_LevelSetObjectBC(size_t field,
                       LevelSetObject* lso,
                       size_t abuse_field,
                       OP op);

  virtual void operator ()();

};

template <typename OP>
CPU_LevelSetObjectBC<OP>::CPU_LevelSetObjectBC(size_t field,
                                               LevelSetObject* lso,
                                               size_t abuse_field,
                                               OP op)
  : LevelSetObjectBC(field, lso, abuse_field)
{
  //  m_Field = field;
  //  m_LevelSetObject = lso;
  //  m_AbuseField = abuse_field;
  m_Op = op;
}

template <typename OP>
void CPU_LevelSetObjectBC<OP>::operator()()
{

  // Potential recursion: Interpolate sets may contain cells in m_InnerCellsLayers.
  // To prevent recursion, the abuse_field is used as intermediate data storage.

  real relax = 0.5;

  /** @todo a better relaxation method would help: Reduce/omit unphysical influence of
    *       relax < 1. , if corrections are small, but keep it to react on inpulsive
    *       starts or similar. Found approx. 0.8 to be a limit for impulsive starts. */

  // Inner Cells
  //.. 1st loop: acquire data, avoid recursion
  for (size_t ll_c = 0; ll_c < m_NumInnerLayerCells; ll_c++) {
    m_Op.access(m_InnerCellsLayers[ll_c], m_Field, m_AbuseField);
  }
  //.. 2nd loop: set values
  for (size_t ll_c = 0; ll_c < m_NumInnerLayerCells; ll_c++) {
    m_Op.operateInner(m_InnerCellsLayers[ll_c], m_Field, m_AbuseField, relax);
  }

  // Outer Cells
  //.. 1st loop: acquire data, avoid recursion
  for (size_t ll_c = 0; ll_c < m_NumOuterLayerCells; ll_c++) {
    m_Op.access(m_OuterCellsLayers[ll_c], m_Field, m_AbuseField);
  }
  //.. 2nd loop: set values
  for (size_t ll_c = 0; ll_c < m_NumOuterLayerCells; ll_c++) {
    m_Op.operateOuter(m_OuterCellsLayers[ll_c], m_Field, m_AbuseField, relax);
  }
}

#endif // CPU_LEVELSETOBJECTBC_H
