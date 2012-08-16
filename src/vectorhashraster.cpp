
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

