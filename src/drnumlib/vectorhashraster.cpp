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

template <class T>
VectorHashRaster<T>::VectorHashRaster ()
  : CartesianRaster()
{
  // nothing ?
}


template <class T>
void VectorHashRaster<T>::setUp (real m_x_min, real m_y_min, real m_z_min,
                                 real m_x_max, real m_y_max, real m_z_max,
                                 size_t num_i, size_t num_j, size_t num_k)
{
  setupAligned(m_x_min, m_y_min, m_z_min,
               m_x_max, m_y_max, m_z_max);
  resize(num_i, num_j, num_k);
  m_data.resize(num_i * num_j * num_k);
  // Clear contents of 2nd dim vectors
  for (size_t l = 0; l < m_data.size(); l++) {
    m_data[l].clear();
  }
}

