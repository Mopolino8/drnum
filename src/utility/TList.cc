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


//#include <typeinfo>

//bisher nicht bemueht macht evtl keinen sinn!!!!
template <class T>
TList<T>::TList (T a_default_value) : List()
{
  value = NULL;
  default_value = a_default_value;
  found_last_time = false;
}

template <class T>
TList<T>::TList (size_t a_max_num_entries, size_t a_delta_entries, T a_default_value) 
  : List (a_max_num_entries, a_delta_entries)
{
  size_t mne = MaxNumEntries();
  value = new T [mne];
  default_value = a_default_value;
  found_last_time = false;
}

template <class T>
TList<T>::TList (List *a_master, T a_default_value) : List (a_master)
{
  //cout<<"TList<T>::TList (List *a_master, T a_default_value) : List (a_master)"<<endl;
  size_t mne = MaxNumEntries();
  value = new T [mne];
  default_value = a_default_value;
  InitAllEntries(default_value);
  found_last_time = false;
}

// gehts ????
template <class T>
TList<T>::TList(const TList<T> &other, const T a_default_value) : List()
{
  value = NULL;
  default_value = a_default_value;
  found_last_time = false;
  operator=(other);
}

template <class T>
void TList<T>::operator=(const TList<T> &other)
{
  List::operator=(other);
  for(size_t i = begin_idx(); i < end_idx(); i = next_idx(i)) {
    value[i] = other[i];
  };
}

template <class T>
TList<T>::~TList () {
  if (value) {
//      for (size_t i = 0; i < MaxNumEntries(); i++) {
//	  delete value[i];
//      };
      delete [] value;
  }
}  

//
//.. InitAllEntries
//
template <class T>
void TList<T>::InitAllEntries (T a_value) {
  size_t i;
  for (i = 0; i < MaxNumEntries(); i++) value[i] = a_value;
}

//
//.. Extend
//
template <class T>
void TList<T>::Extend (size_t delta) {
  T *new_value;
  size_t i_entry;
  new_value = new T [MaxNumEntries() + delta];
  for (i_entry = 0; i_entry < MaxNumEntries(); i_entry++) new_value[i_entry] = value[i_entry];
  delete [] value;
  value = new_value;
}

//
//.. FindItem
//
template <class T>
size_t TList<T>::FindItem (T item)
{
  if (found_last_time) {
    if (IsActive(last_found)) {
      if (value[last_found] == item) return last_found;
    };
  };
  for(size_t i = begin_idx(); i < end_idx(); i = next_idx(i)) {
    if (value[i] == item) {
      found_last_time = true;
      last_found = i;
      return i;
    };
  };
  string msg = "TList::FindItem\n";
  msg += typeid(*this).name();
  throw NotFound_error(msg);
}

//
//.. IsIn
//
template <class T>
bool TList<T>::Has(T item)
{
  if (found_last_time) {
    if (IsActive(last_found)) {
      if (value[last_found] == item) return true;
    };
  };
  for(size_t i = begin_idx(); i < end_idx(); i = next_idx(i)) {
    if (value[i] == item) {
      found_last_time = true;
      last_found = i;
      return true;
    };
  };
  return false;
}
/*
// keine Ahnung wozu benoetigt ATTENTION
template <class T>
ostream& operator<<(ostream &s, const TList<T> &tlist)
{
  s << '[';
  FORALL(i, tlist.) {
    s << "(" << i << "," << tlist[i] << ")";
    if (i != tlist.last_idx()) s << ',';
  };
  s << ']';
  return s;
}
*/ 

/*             
// geaendert am 09.10.06 old version was wrong
template <class T>
void TList<T>::SortUp()
{
  // geht auch in die Hose wegen der Definition von last_idx() 
  size_t i_start = begin_idx();
  if(begin_idx() != end_idx())
  for (size_t i = i_start; i < last_idx(); i = next_idx(i)) 
    {
      for (size_t j = i_start; j < last_idx(); j = next_idx(j)) 
	{
	  if (At(i) > At(next_idx(i))) 
	    {
	      Swap(i, next_idx(i));
	    };
	};
    };
}
// geaendert am 09.10.06  old version was wrong
template <class T>
void TList<T>::SortDown()
{
  // geht auch in die Hose wegen der Definition von last_idx() 
  size_t i_start = begin_idx();
  if(begin_idx() != end_idx())
  for (size_t j = i_start; j < last_idx(); j = next_idx(j)) 
    {
      for (size_t i = i_start; i < last_idx(); i = next_idx(i)) 
	{
	  if (At(i) < At(next_idx(i))) 	Swap(i, next_idx(i));
	};
    };
}
*/



// nochmal geaendert am 09.10.06 old version was wrong
template <class T>
void TList<T>::SortUp()
{
    size_t i_start = begin_idx();
    bool proceed = true;
    while(proceed){
	proceed = false;
	for(size_t i = i_start; i < end_idx(); i = next_idx(i)){
	    size_t i_next = next_idx(i);
	    if(i_next != end_idx()) {  // Schleife leider nicht abzaehlbar
		if((At(i) > At(i_next))) {
		    Swap(i, i_next);
		    proceed = true;
		};
	    };
	};
    };
};
// nochmal geaendert am 09.10.06 old version was wrong
template <class T>
void TList<T>::SortDown()
{
    size_t i_start = begin_idx();
    bool proceed = true;
    while(proceed){
	proceed = false;
	for(size_t i = i_start; i < end_idx(); i = next_idx(i)){
	    size_t i_next = next_idx(i);
	    if(i_next != end_idx()) {  // Schleife leider nicht abzaehlbar
		if((At(i) < At(i_next))) {
		    Swap(i, i_next);
		    proceed = true;
		};
	    };
	};
    };
};

template <class T>
T TList<T>::MinValue()
{
  T min_value = At(0);
  FORALL(i, this->) {
    if(At(i) < min_value) {
      min_value = At(i);
    };
  };
  return min_value;
}

template <class T>
T TList<T>::MaxValue()
{
  T max_value = At(0);
  FORALL(i, this->) {
    if(max_value < At(i)) {
      max_value = At(i);
    };
  };
  return max_value;
}
