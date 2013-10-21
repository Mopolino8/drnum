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
#include "structuredhexraster.h"

StructuredHexRaster::StructuredHexRaster()
{
  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
  m_NumL = m_NumI * m_NumJ * m_NumK; // = 1
  //m_CellLink = new List(1,1);
  //m_CellLink = NULL;
  m_Eps = 1.e-6;  /// @todo need a better eps-handling.
}


void StructuredHexRaster::resize(const size_t& num_i, const size_t& num_j, const size_t& num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  m_NumL = m_NumI * m_NumJ * m_NumK; // = 1
}


