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
#if !defined(SPARSEWEIGHTEDSET_HH)
#define SPARSEWEIGHTEDSET_HH

/**
  * Sparse version of WeightedSet with restricted functions. Ommit stl::vector, to avoid
  * vast memory overhead.
  *
  * Functions:
  * - copy from a WeightedSet<T>
  * - multiply weights with a scalar
  * - calculate data
**/
template<class T> class SparseWeightedSet;

#include <algorithm>
#include <vector>
#include <iostream>

#include "weightedset.h"

using namespace std;

/**
  * Sparse version of SparseWeightedSet with restricted functions. Ommit stl::vector, to avoid
  * vast memory overhead.
  *
  * Functions:
  * - copy from a WeightedSet<T>
  * - multiply weights with a scalar
  * - calculate data
**/
template<class T>
class SparseWeightedSet
{
protected: // attributes
  /// @todo need usefull protection rules

  pair<size_t, T>* v;  // NOTE: pair smaller than single lists
  // size_t* jj;
  // T* ww;

  size_t m_Size;

  size_t m_Capacity;

  // T t_zero;
  // T t_one;

protected: // methods

  void ensureCapacity(const size_t& min_size);  ///< increases capacity to fit at least min_size
  void adaptSize();                             ///< adapt size of v to m_Capacity
  void setup();                                 ///< set t_zero and t_one

public:

  /** Dummy construction
  */
  // Constructors
  SparseWeightedSet();      ///< Dummy construction, empty contents

  // Member functions
  bool reserve(const size_t& size);               ///< reserves space for insertions
  void pushBack(const size_t& j, const T& t);     ///< appends a pair (j,s)
  void clear();                                   ///< Erase all contents
  void fitToSize();                               ///< Reduce capacity to size needed
  void PrintNumSets();                            ///< print to screen
  void MultScalar(const T& s);                    ///< multiply     v = s*v
  T WeightSum() const;                            ///< compute sum of all weights
  size_t getSize() const {return m_Size;}         ///< return number of pairs.
  size_t getCapacity() const {return m_Capacity;} ///< return reserved size
  size_t getIndex(const size_t& ll) const;        ///< return the index of the ll-th pair (index, weight)
  T getWeight(const size_t& ll) const;            ///< return the weight of the ll-th pair (index, weight)
  T computeValue(const T* a) const;               ///< Compute sum(a[j[i]] * t[i]) for all i

  // Operators
  void operator=(const SparseWeightedSet<T>& w);  ///< copy:      v = w
  void operator=(const WeightedSet<T>& a_ws);     ///< copy:      v = w
  void operator*=(const T& s);                    ///< multiply:  v = s*v

};


template<class T>
inline void SparseWeightedSet<T>::ensureCapacity(const size_t& min_size)
{
  if(m_Capacity < min_size) {
    if(m_Capacity == 0) {
      m_Capacity = 1;
    } else {
      m_Capacity *= 2;
    }
  }
  adaptSize();
}

template<class T>
inline void SparseWeightedSet<T>::adaptSize()
{
//  size_t* jj_h = new size_t[m_Capacity];
//  T* ww_h = new T[m_Capacity];
//  for (size_t i = 0; i < m_Size; i++) {
//    jj_h[i] = jj[i];
//    ww_h[i] = ww[i];
//  }
//  delete[] jj;
//  jj = jj_h;
//  delete[] ww;
//  ww = ww_h;
  pair<size_t, T>* v_h = new pair<size_t, T>[m_Capacity];
  for (size_t i = 0; i < m_Size; i++) {
    v_h[i] = v[i];
  }
  delete[] v;
  v = v_h;
}

template<class T>
inline void SparseWeightedSet<T>::setup()
{
  //t_zero = T(0);
  //t_one = T(1);
}

template<class T>
inline SparseWeightedSet<T>::SparseWeightedSet()
{
  // jj = NULL;
  // ww = NULL;
  v = NULL;
  setup();
  clear();
}

template<class T>
inline bool SparseWeightedSet<T>::reserve(const size_t& size)
{
  bool error = m_Capacity > size; // will not shrink
  if(size > m_Capacity) {
    m_Capacity = size;
    adaptSize();
    error = false;
  }
  return error;
}

template<class T>
inline void SparseWeightedSet<T>::pushBack(const size_t& j, const T& s) {
  ensureCapacity(m_Size+1);
  // jj[m_Size] = j;
  // ww[m_Size] = s;
  pair<size_t, T> new_pair;
  new_pair.first = j;
  new_pair.second = s;
  v[m_Size] = new_pair;
  m_Size++;
}

template<class T>
inline void SparseWeightedSet<T>::clear()
{
//  if(jj) {
//    delete[] jj;
//  }
//  if(ww) {
//    delete[] ww;
//  }
  if(v) {
    delete[] v;
  }
  m_Size = 0;
  m_Capacity = 0;
}

template<class T>
inline void SparseWeightedSet<T>::fitToSize()
{
  if(m_Capacity != m_Size) {
    m_Capacity = m_Size;
    adaptSize();
  }
}

template<class T>
inline void SparseWeightedSet<T>::PrintNumSets()
{
  T weight_sum = WeightSum();
  for(size_t i=0; i<m_Size; i++) {
     cout << "i = " << v[i].first << " ; weight = " << v[i].second  << " rel_weight = " << (v[i].second/weight_sum) << endl;
//    cout << "i = " << jj[i] << " ; weight = " << ww[i] << " rel_weight = " << (ww[i]/weight_sum) << endl;
  }
  cout << "sum of weights = " << weight_sum << endl;
}

template<class T>
inline void SparseWeightedSet<T>::MultScalar(const T& scalar)
{
  for (size_t i = 0; i < m_Size; i++) {
    // ww[i] *= scalar;
    v[i].second *= scalar;
  }
}

template<class T>
inline T SparseWeightedSet<T>::WeightSum() const
{
  T weight_sum = T(0);
  //T weight_sum = t_zero;
  for (size_t i = 0; i < m_Size; i++) {
    weight_sum += v[i].second;
    // weight_sum += ww[i];
  }
  return weight_sum;
}

template<class T>
inline size_t SparseWeightedSet<T>::getIndex(const size_t& ll) const
{
  return v[ll].first;
  // return jj[ll];
}

template<class T>
inline T SparseWeightedSet<T>::getWeight(const size_t& ll) const
{
  return v[ll].second;
  // return ww[ll];
}

template<class T>
inline T SparseWeightedSet<T>::computeValue(const T* a) const
{
  T ret_val = T(0);  /// @todo Be sure, this will not cause a performance issue!!!
  //T weight_sum = t_zero;
  for (size_t i = 0; i < m_Size; i++) {
    ret_val += a[v[i].first] * v[i].second;
    // ret_val += a[jj[i]] * ww[i];
  }
  return ret_val;
}

template<class T>
inline void SparseWeightedSet<T>::operator=(const SparseWeightedSet<T>& a_ws)
{
  clear();
  reserve(a_ws.getSize());  // also sets capacity
  m_Size = m_Capacity;
  for (size_t i = 0; i < m_Size; i++) {
    v[i] = a_ws.v[i];
    // jj[i] = a_ws.jj[i];
    // ww[i] = a_ws.ww[i];
  }
}

template<class T>
inline void SparseWeightedSet<T>::operator=(const WeightedSet<T>& a_ws)
{
  clear();
  reserve(a_ws.getSize());
  m_Size = m_Capacity;
  for (size_t i = 0; i < m_Size; i++) {
    v[i] = a_ws.v[i];
    // jj[i] = a_ws.v[i].first;
    // ww[i] = a_ws.v[i].second;
  }
}

template<class T>
inline void SparseWeightedSet<T>::operator*=(const T& scalar) {
  MultScalar(scalar);
}

#endif
