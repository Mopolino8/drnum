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
#if !defined(USPARSEWEIGHTEDSET_HH)
#define USPARSEWEIGHTEDSET_HH

/**
  * Sparse version of WeightedSet with restricted functions. Ommit stl::vector, to avoid
  * vast memory overhead.
  *
  * Functions:
  * - copy from a WeightedSet<T>
  * - multiply weights with a scalar
  * - calculate data
**/
template<class T> class USparseWeightedSet;

#include <algorithm>
#include <vector>
#include <iostream>

#include "weightedset.h"

using namespace std;

/**
  * Ultra sparse version of WeightedSet with restricted functions. Ommit stl::vector, to avoid
  * vast memory overhead. Also discard efficient growth mechanism to save capacity variable.
  *
  * Functions:
  * - copy from a WeightedSet<T>
  * - multiply weights with a scalar
  * - calculate data
**/
template<class T>
class USparseWeightedSet
{
protected: // attributes
  /// @todo need usefull protection rules

  pair<size_t, T>* v;  // NOTE: pair smaller than single lists
  // size_t* jj;
  // T* ww;

  size_t m_Size;

  // size_t m_Capacity;

protected: // methods

  void setCapacity(const size_t& size);  ///< sets capacity to size

public:

  /** Dummy construction
  */
  // Constructors
  USparseWeightedSet();      ///< Dummy construction, empty contents

  // Member functions
  void pushBack(const size_t& j, const T& t);     ///< appends a pair (j,s)
  void clear();                                   ///< Erase all contents
  void PrintNumSets();                            ///< print to screen
  void MultScalar(const T& s);                    ///< multiply     v = s*v
  T WeightSum() const;                            ///< compute sum of all weights
  size_t getSize() const {return m_Size;}         ///< return number of pairs.
  size_t getIndex(const size_t& ll) const;        ///< return the index of the ll-th pair (index, weight)
  T getWeight(const size_t& ll) const;            ///< return the weight of the ll-th pair (index, weight)
  T computeValue(const T* a) const;               ///< Compute sum(a[j[i]] * t[i]) for all i

  // Operators
  void operator=(const USparseWeightedSet<T>& w);  ///< copy:      v = w
  void operator=(const WeightedSet<T>& a_ws);     ///< copy:      v = w
  void operator*=(const T& s);                    ///< multiply:  v = s*v

};


template<class T>
inline void USparseWeightedSet<T>::setCapacity(const size_t& size)
{
  /// @todo this is critical as it leaves m_Size unequal to actual array allocation
  pair<size_t, T>* v_h = new pair<size_t, T>[size];
  m_Size = min(size, m_Size);
  for (size_t i = 0; i < m_Size; i++) {
    v_h[i] = v[i];
  }
  if(v) {
    delete[] v;
  }
  v = v_h;
}

template<class T>
inline USparseWeightedSet<T>::USparseWeightedSet()
{
  // jj = NULL;
  // ww = NULL;
  v = NULL;
  clear();
}

template<class T>
inline void USparseWeightedSet<T>::pushBack(const size_t& j, const T& s) {
  // NOTE: will disrespect any previous capacity setting (setCapacity)
  setCapacity(m_Size+1);
  pair<size_t, T> new_pair;
  new_pair.first = j;
  new_pair.second = s;
  v[m_Size] = new_pair;
  m_Size++;
}

template<class T>
inline void USparseWeightedSet<T>::clear()
{
  if(v) {
    delete[] v;
  }
  m_Size = 0;
}

template<class T>
inline void USparseWeightedSet<T>::PrintNumSets()
{
  T weight_sum = WeightSum();
  for(size_t i=0; i<m_Size; i++) {
     cout << "i = " << v[i].first << " ; weight = " << v[i].second  << " rel_weight = " << (v[i].second/weight_sum) << endl;
//    cout << "i = " << jj[i] << " ; weight = " << ww[i] << " rel_weight = " << (ww[i]/weight_sum) << endl;
  }
  cout << "sum of weights = " << weight_sum << endl;
}

template<class T>
inline void USparseWeightedSet<T>::MultScalar(const T& scalar)
{
  for (size_t i = 0; i < m_Size; i++) {
    // ww[i] *= scalar;
    v[i].second *= scalar;
  }
}

template<class T>
inline T USparseWeightedSet<T>::WeightSum() const
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
inline size_t USparseWeightedSet<T>::getIndex(const size_t& ll) const
{
  return v[ll].first;
  // return jj[ll];
}

template<class T>
inline T USparseWeightedSet<T>::getWeight(const size_t& ll) const
{
  return v[ll].second;
  // return ww[ll];
}

template<class T>
inline T USparseWeightedSet<T>::computeValue(const T* a) const
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
inline void USparseWeightedSet<T>::operator=(const USparseWeightedSet<T>& a_ws)
{
  clear();
  setCapacity(a_ws.getSize());
  m_Size = a_ws.getSize();
  for (size_t i = 0; i < m_Size; i++) {
    v[i] = a_ws.v[i];
    // jj[i] = a_ws.jj[i];
    // ww[i] = a_ws.ww[i];
  }
}

template<class T>
inline void USparseWeightedSet<T>::operator=(const WeightedSet<T>& a_ws)
{
  clear();
  setCapacity(a_ws.getSize());
  m_Size = a_ws.getSize();
  for (size_t i = 0; i < m_Size; i++) {
    v[i] = a_ws.v[i];
    // jj[i] = a_ws.jj[i];
    // ww[i] = a_ws.ww[i];
  }
}

template<class T>
inline void USparseWeightedSet<T>::operator*=(const T& scalar) {
  MultScalar(scalar);
}

#endif
