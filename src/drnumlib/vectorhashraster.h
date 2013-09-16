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
#ifndef VECTORHASHRASTER_HH
#define VECTORHASHRASTER_HH

#include <cstddef>
//#include <namespace_mouse.hh>

template <class T> class VectorHashRaster;

#include "cartesianraster.h"
#include <vector>


/** Class to hold a VectorHashRaster on a CartesianRaster
 *
 *  Hold arbitrary ammounts of data types T on each cell of the CartesianRaster
 *
 */
template <class T>
class VectorHashRaster : public CartesianRaster
{
  protected:
  /**
    * vector of vectors to hold data.
    */
  vector<vector<T> > m_data;


  public:
  /** constructor
    */
  VectorHashRaster();


  /**
    * Set up format and resolution.
    */
  void setUp(real m_x_min, real m_y_min, real m_z_min,
             real m_x_max, real m_y_max, real m_z_max,
             size_t num_i, size_t num_j, size_t num_k);


  /**
    * Insert an element at position (i,j,k)
    * @param i i-index
    * @param j j-index
    * @param k k-index
    * @param value value to insert
    */
  void insert(size_t i, size_t j, size_t k, T value);

  /**
    * Insert an element at position l = f(i,j,k)
    * @param l 1D-index in raster
    * @param value value to insert
    */
  void insert(size_t l, T value);

  /**
    * Save insertion of an element at position (i,j,k).
    * Will not insert, if out of bounds.
    * @param i i-index
    * @param j j-index
    * @param k k-index
    * @param value value to insert
    * @return bool indicating insertion was in valid i,j,k bounds
    */
  bool save_insert(size_t i, size_t j, size_t k, T value);

  /**
    * Save insertion of an element at position (i,j,k).
    * Will not insert, if out of bounds.
    * @param l 1D-index in raster
    * @param value value to insert
    * @return bool indicating insertion was in valid i,j,k bounds
    */
  bool save_insert(size_t l, T value);

  /**
    * Make stored entries unique by sorting and erasing.
    * NOTE: Original order will lost due to sorting.
    * @param i i-index
    * @param j j-index
    * @param k k-index
    */
  void unify(const size_t& i, const size_t& j, const size_t& k);

  /**
    * Make stored entries unique by sorting and erasing.
    * NOTE: Original order will lost due to sorting.
    * @param l 1D-index in raster
    */
  void unify(const size_t& l);

  /// @todo No save version of "unify" available.

  /**
    * Get the number of entries in the second dimension of t_insect.
    * @param i i-index
    * @param j j-index
    * @param k k-index
    * @return the number of items stored for l
    */
  size_t getNumItems(const size_t& i, const size_t& j, const size_t& k) const;

  /**
    * Get the number of entries in the second dimension of t_insect.
    * @param l the index of the cell in the raster
    * @return the number of items stored for l
    */
  size_t getNumItems(const size_t& l) const;

  /**
    * Data access method (read/write)
    * @param l 1D-index in raster
    * @param l_i the l_i-th entry for a raster node "l"
    * @return a reference to the entry
    */
  T& at(const size_t& l, const size_t& l_i);
  //{return m_data[l][l_i];};


  /**
    * Data access method (read/write)
    * @param i i-index
    * @param j j-index
    * @param k k-index
    * @param l_i the l_i-th entry for a raster node "l"
    * @return a reference to the entry
    */
  T& at(const size_t& i, const size_t& j, const size_t& k, const size_t& l_i);

  /// @todo No save version of "at" available.

  /** destructor */
  virtual ~VectorHashRaster() {}  /// @todo destructor missing

};


template <class T>
inline void VectorHashRaster<T>::insert(size_t i, size_t j, size_t k, T value)
{
  // insert "value" at given position
  size_t l = index(i, j, k);
  m_data[l].push_back(value);
  unify(l);
}


template <class T>
inline void VectorHashRaster<T>::insert(size_t l, T value)
{
  m_data[l].push_back(value);
  unify(l);
}


template <class T>
inline bool VectorHashRaster<T>::save_insert(size_t i, size_t j, size_t k, T value)
{
  bool error = false;
  size_t l = save_index(i, j, k,
                        error);
  if (!error) {
    m_data[l].push_back(value);
    unify(l);
  }
  return !error;
}


template <class T>
inline bool VectorHashRaster<T>::save_insert(size_t l, T value)
{
  bool error = l > sizeL();
  if (!error) {
    m_data[l].push_back(value);
    unify(l);
  }
  return !error;
}


template <class T>
inline void VectorHashRaster<T>::unify(const size_t& i, const size_t& j, const size_t& k)
{
  size_t l = index(i, j, k);
  sort(m_data[l].begin(), m_data[l].end());
  // Remove duplicates
  typename vector<T>::iterator it;
  it = unique(m_data[l].begin(), m_data[l].end());
  m_data[l].resize(it - m_data[l].begin());
}


template <class T>
inline void VectorHashRaster<T>::unify(const size_t& l)
{
  sort(m_data[l].begin(), m_data[l].end());
  // Remove duplicates
  typename vector<T>::iterator it;
  it = unique(m_data[l].begin(), m_data[l].end());
  m_data[l].resize(it - m_data[l].begin());
}


template <class T>
inline size_t VectorHashRaster<T>::getNumItems(const size_t& i, const size_t& j, const size_t& k) const {
  size_t l = index(i, j, k);
  return m_data[l].size();
}


template <class T>
inline size_t VectorHashRaster<T>::getNumItems(const size_t& l) const {
  return m_data[l].size();
}


template <class T>
inline T& VectorHashRaster<T>::at(const size_t& i, const size_t& j, const size_t& k, const size_t& l_i) {
  size_t l = index(i, j, k);
  return m_data[l][l_i];
}


template <class T>
inline T& VectorHashRaster<T>::at(const size_t& l, const size_t& l_i) {
  return m_data[l][l_i];
}



//{return m_data[l][l_i];};

//template <class T>
//inline T VectorHashRaster<T>::at(const size_t& l, const size_t& l_i)
//{
//  return m_data[l][l_ị];
//}


#include "vectorhashraster.cpp"

#endif  //VectorHashRaster_HH
