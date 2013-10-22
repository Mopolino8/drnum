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

#include "drnum.h"

protected: // attributes

size_t m_Field;                   /// the variable field to work on
size_t m_AbuseField;              /// abused field to avoid recursion

size_t m_NumInnerLayerCells;      /// total number of cells in inner layers
size_t m_NumOuterLayerCells;      /// total number of cells in outer layers

/** Storage of levelset layer data for inner cells
  * Data alignment (all with actual length, no padding)
  * 1st dim: patch id (direct; 0, 1, 2, ...)
  * 2nd dim: layer:
  *       0 :  has at least one face neighbour with negative G-value (inside)
  *       1 .. m_NumOuterLayers : furter layers farther outside
  * 3rd dim: cell indices in layer */
LSLayerDataExtrapol* m_InnerCellsLayers;

/** Storage of levelset layer data for outer cells
  * Data alignment (all with actual length, no padding)
  * 1st dim: patch id (direct; 0, 1, 2, ...)
  * 2nd dim: layer:
  *       0 :  has at least one face neighbour with negative G-value (inside)
  *       1 .. m_NumOuterLayers : furter layers farther outside
  * 3rd dim: cell indices in layer */
LSLayerDataExtrapol* m_OuterCellsLayers;
