// !!
// This is a part of MOUSE, a library for PDE's on unstructured grids
// Copyright (C) 1999 Oliver Gloth <oliver@vug.uni-duisburg.de>
// Institut fuer Verbrennung und Gasdynamik (Universitaet Duisburg)
// institute for combustion and gas dynamics (university of Duisburg)
// Thursday, 1 April, 1999 Duisburg, Germany
//
// please see http://www.vug.uni-duisburg.de/MOUSE for more information
// please send any questions or suggestions to mouse@www.vug.uni-duisburg.de
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// !!

#if !defined(TLIST_H)
#define TLIST_H

#include <cstddef>
//#include <typeinfo>

#include <namespace_mouse.hh>
BEGIN_MOUSE
template <class T> class TList;
template <class T, size_t MNE, size_t DE> class TTList;
END_MOUSE

#include <List.hh>
#include <typeinfo>


BEGIN_MOUSE
/**
 * a self-allocating one-dimensional array.
 */
template <class T>
class TList : public List 
{
private:
  /// a flag which is true if Has or FindItem has already been successfully called
  bool found_last_time;

  /// the index of the last item found
  size_t last_found;

protected:
  /** array of the stored values */
  T *value;
  
  /** extends the array value.
      @param delta the increment */ 
  virtual void Extend (size_t delta);

  /** initialize an entry with the default value 
      @param i the entry to initialize */ 
  virtual void InitEntry (size_t i) {
    value[i] = default_value;
  };
  
  virtual void CopyEntry (size_t source_point, size_t dest_point) {
    List::CopyEntry(source_point, dest_point);
    value[dest_point] = value[source_point];
  };

public:
  /// the default value for new entries
  T default_value;
    
  /**
   * default constructor.
   * @param a_default_value for new entries 
   */
  TList(T a_default_value = T());

  /** 
   * @param a_max_num_entries the initial maximal number of entries
   * @param a_delta_entries the increment for expansion
   * @param a_default_value for new entries 
   */
  TList(size_t a_max_num_entries , size_t a_delta_entries , T a_default_value = T());
    
  /**
   * @param a_master a master for this List.
   * @param a_default_value the initial value for new entries 
   */
  TList(List *a_master, T a_default_value = T());

  TList(const TList<T> &other, const T a_default_value = T());
  
  virtual ~TList ();
  
  /**
   * Initialize all active entries with a value.
   * @param v the value for the initilaziation 
   */
  void InitAllEntries (T v);

  /**
   * Allocate a new entry.
   * @return a reference to the new entry 
   */
  T& NewEntry () { 
    size_t i = AddEntry ();  
    return value[i];  
  };

  /**
   * Insert an entry at given address i
   * @param i the address to insert in
   * @return a reference to the new entry 
   */
  T& InsertEntry (const size_t& i) {
//      if(i >= NumEntries()) {
//	  size_t growth_needed = i - NumEntries() + 1;
//	  Alloc(growth_needed);
//      };
      while(i >= NumEntries()) {  // ATTENTION: better use List::Alloc(growth_needed) later
	  AddEntry();
      };
      return value[i];
  }

  /** The Value of an entry.
   *  @param i the index of the entry
   *  @return a reference to the value of this entry 
   */
//  T& At(size_t i) {return value[i];}
    T& At (const size_t& i) { 
#ifdef MDEBUG
	if (i >= end_idx()) {  // note, that negative input would be permutated in large size_t too.
	    cerr << "TList::T& At (size_t i): index " << i << " out of bounds" << endl;
	    throw InvalidIndex_error(i);
	};
	if (!IsActive(i)) {
	    cerr << "TList::T& At (size_t i): the entry number " << i << " is inactive" << endl;
	    throw InvalidIndex_error(i);
	};
#endif
	return value[i]; 
    };
    
    /** The Value of an entry.
     *  @param i the index of the entry
     *  @return the value of this entry
     */
    T const& At (size_t i) const { 
#ifdef MDEBUG
	if (i >= end_idx()) {
	    cerr << "TList::T const& At (size_t i) const: index " << i << " out of bounds" << endl;
	    throw InvalidIndex_error(i);
	};
	if (!IsActive(i)) {
	    cerr << "TList::T const& At (size_t i) const: the entry number " << i << " is inactive" << endl;
	    throw InvalidIndex_error(i);
	};
#endif
	return value[i]; 
    };
    
    /** Reference-argument version to access value. Unfiortunately, no performance
     *  impact was found on gcc and intel comilers.
     *  @param i the index of the entry
     *  @return the value of this entry
     */
    T& AtR (size_t& i) { 
#ifdef MDEBUG
	if (i >= end_idx()) {
	    cerr << "TList::T& At (size_t i): index " << i << " out of bounds" << endl;
	    throw InvalidIndex_error(i);
	};
	if (!IsActive(i)) {
	    cerr << "TList::T& At (size_t i): the entry number " << i << " is inactive" << endl;
	    throw InvalidIndex_error(i);
	};
#endif
	return value[i]; 
    };

  /** The Value of an entry.
      @param i the index of the entry
      @return a reference to the value of this entry */
  T& operator[] (size_t i)  { return At(i); }

  /** The Value of an entry.
      @param i the index of the entry
      @return the value of this entry */
  const T operator[] (size_t i) const { return At(i); }

  /** Find an entry.
    @param v the value to search for
    @return the index of the entry */
  size_t FindItem (T v);

  /** Check if an entry is in the list */
  bool Has(T v);

  /** Get a pointer to the data field */
  T* ValuePointer() { return value; };

  /** exchange two entries
   *  @param i1 the first entry
   *  @param i2 the second entry
   */ 
  virtual void Swap(size_t i1, size_t i2) {
    T h = At(i1);
    At(i1) = At(i2);
    At(i2) = h;
  };

  /** */ void operator=(const TList<T> &other);

  /** return the first active entry
   *  @return a reference to the first active entry.
   */
  // WARNING auskommentiert weil nicht deklariert
  //   T& FirstEntry() { return At(first_idx()); };

  /** return the first active entry
   *  @return a reference to the last active entry.
   */
  T& LastEntry() { return At(last_idx()); };

  /** Output operator.
   */ 

  // ATTENTION  superfluid ?????
   template <class aT> 
   friend ostream& operator<<(ostream &s, const TList<aT> &tlist);

  

  // =====================================
  // STL style iterator and related things
  // =====================================

  typedef T value_type;

  class iterator {

    TList<T> *tlist;
    size_t i_current;

  public:
    iterator(TList<T> *a_tlist, size_t i) : tlist(a_tlist), i_current(i) {};
    iterator() : tlist(NULL), i_current(0) {};
    bool operator==(const iterator &iter) const;
    iterator& operator++();
    iterator operator++(int);
    T& operator*();
  };

  class const_iterator {

    const TList<T> *tlist;
    size_t i_current;

  public:
    const_iterator(const TList<T> *a_tlist, size_t i) : tlist(a_tlist), i_current(i) {};
    const_iterator() : tlist(NULL), i_current(0) {};
    bool operator==(const const_iterator &iter) const;
    const_iterator& operator++();
    const_iterator operator++(int);
    T operator*();
  };

  iterator begin() { return iterator(this, begin_idx()); };
  iterator end() { return iterator(this, end_idx()); };
  const_iterator begin() const { return const_iterator(this, begin_idx()); };
  const_iterator end() const { return const_iterator(this, end_idx()); };

  /**
   * sort the TList.
   * This only works for types <it>T</it> which have
   * the operator '>' and '='. The smallest item will be
   * the first after this sorting.
   */
  void SortUp();

  /**
   * sort the TList.
   * This only works for types <it>T</it> which have
   * the operator '<' and '='. The smallest item will be
   * the last after this sorting.
   */
  void SortDown();

  /**
   * Find minimum value of all entries.
   * This only works for types <it>T</it> which have
   * the operator '<' and '='.
   */
  T MinValue();

  /**
   * Find maximum value of all entries.
   * This only works for types <it>T</it> which have
   * the operator '<' and '='.
   */
  T MaxValue();

};


// TList<T>::iterator
// ==================
template<class T>
inline typename TList<T>::iterator& TList<T>::iterator::operator++()
{ 
  i_current = tlist->next_idx(i_current);
  return *this;
}

template<class T>
inline typename TList<T>::iterator TList<T>::iterator::operator++(int)
{
  size_t old_i_current;
  i_current = tlist->next_idx(i_current);
  return iterator(tlist, old_i_current); 
}

template<class T>
inline bool TList<T>::iterator::operator==(const iterator &iter) const
{
  return ((iter.tlist == tlist) && (iter.i_current == i_current));
}

template<class T>
inline T& TList<T>::iterator::operator*()
{
  return (*tlist)[i_current];
}

// TList<T>::const_iterator
// ========================
template<class T>
inline typename TList<T>::const_iterator& TList<T>::const_iterator::operator++()
{ 
  i_current = tlist->next_idx(i_current);
  return const_iterator(tlist, i_current); 
}

template<class T>
inline typename TList<T>::const_iterator TList<T>::const_iterator::operator++(int)
{
  size_t old_i_current;
  i_current = tlist->next_idx(i_current);
  return const_iterator(tlist, old_i_current); 
}

template<class T>
inline bool TList<T>::const_iterator::operator==(const const_iterator &iter) const
{
  return ((iter->tlist == tlist) && (iter.i_current == i_current));
}

template<class T>
inline T TList<T>::const_iterator::operator*()
{
  return (*tlist)[i_current];
}


template<class T, size_t MNE, size_t DE>
class TTList : public TList<T>
{
public:
  TTList() : TList<T>(MNE, DE) {};
};

#include "TList.cc"
END_MOUSE

#endif









