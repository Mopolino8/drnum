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
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <class T>
TInsectHashRaster<T>::TInsectHashRaster (real m_x_min, real m_y_min, real m_z_min,
                                         real m_x_max, real m_y_max, real m_z_max,
                                         size_t num_i, size_t num_j, size_t num_k,
                                         size_t content_size, size_t content_increment)
    : CartesianRaster()
{
  setupAligned(m_x_min, m_y_min, m_z_min,
               m_x_max, m_y_max, m_z_max);
  resize(num_i, num_j, num_k);
  //
  // Build a TInsectionList<t> linked to List* m_CellLink. Main dimension
  // of t_insect gets size equivalent to m_CellLink, this is the number of "cells" in
  // CartesianRaster (inherited) by "this".
  // Size an increment of content lists given by argument or computed.
  if(content_size == 0) { // this is the default
      //.. Set to same size as main list. Average one per cell in CartesianRaster.
      content_size = m_CellLink->NumEntries();
  };
  if(content_increment == 0) { // this is the default
      content_increment = content_size;
  };
  m_tinsect = new TInsectionList<T>(m_CellLink, content_size, content_increment);
  initInsect();
}

