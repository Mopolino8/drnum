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

#if !defined(List_HH)
#define List_HH

#include "namespace_mouse.hh"

BEGIN_MOUSE


class List;


END_MOUSE

#include <cstddef>
//#include <pthread.h>

#include "mouse_types.hh"
#include "MError.hh"

BEGIN_MOUSE

/**
 * This exception will be thrown if an item cannot be found in a List,
 * or one of it's derived classes.
 */
struct NotFound_error : public MError
{
  int code;
  NotFound_error(string msg) : MError(msg) { code = 0; };
  NotFound_error(string msg, int a_code)  : MError(msg) { code = a_code; };
};

/// This exception will be thrown if 'mdebug' is enabled and an index is invalid.
struct InvalidIndex_error : public MError
{
  size_t i;
  InvalidIndex_error(size_t an_i)  : MError("InvalidIndex_error")
  {
    i = an_i;
  };
};

/**
 * This is a dynamic list. It provides only the necessary handling
 * routines and control data, but no payload. Other Lists can be linked to ensure
 * a consistent length of different Lists. They will expand and shrink simultaneously. 
 * Please note that List and its
 * derived classes are not thread-safe in all cases.
 * It is, however, planned to make them thread-save.
 */
class List 
{
public:

  typedef long unsigned int state_t;

private:
  
  /// the first block where new entries might be.
  size_t first_block;

  /// one after the last active_entry.
  size_t after_last_entry;

  /// the dynamical increment.
  size_t delta_entries;

  /// the maximal number of entries.
  size_t max_num_entries;

  /// array to store the entry mapping after garbage collection.
  size_t *offset;

  /// the maximal offset value found.
  size_t max_offset;

  /// indicates if a certain entry is active or not.
  bool *active;

  /// number of client lists if this is a master.
  unsigned int num_clients;

  /// the master List, which is managing the dynamics. 
  List *master;

  /// a field to store pointers to the clients.
  List **client;

  /// flag to determine if the list is clean (no holes).
  bool is_clean;

  /// field for the block sizes.
  size_t *block_size;

  /// update the block_size field to represent current state.
  void UpdateBlockSize();

  /// switch to determine if the block_size field is updated.
  bool block_size_updated;

  /// an internal counter representing list states @see State
  state_t state_counter;

  /**
   * a mutex used if any maintenence function is called.
   * This will only happen in the master List.
   */
  //pthread_mutex_t mutex;

protected:

  /** 
   * initialize the List.
   * @param mne the initial maximum number of entries
   * @param delta a delta entries to extend or shrink the list
   */
  void InitList (size_t mne, size_t delta);

  /**
   * extend or shrink the list.
   * @param delta the number of new entries, negative values will shrink the list
   */ 
  void ExtendList(size_t delta);

  /**
   * add a new client if this is a master list.
   * @param client_to_add a pointer to the new client list
   */ 
  void AddClient (List *client_to_add);
  
  /** 
   * delete a client from this master.
   * @param client_to_del a pointer to the client,
   * which shall be deleted
   */
  void DelClient (List *client_to_del);

  /// flag to determine if the list has been initialized. 
  bool initialized;


  bool IsInitialized () const {return initialized;}

  /**
   * initialization of an entry.
   * This method does nothing for the base-class and can be overwritten in child-classes.
   * @param i the index of the entry 
   */
  virtual void InitEntry (size_t) {};

  /**
   * extend the list.
   * This method has to be overwritten in child-classes to ensure a proper 
   * working of the cleaning mechanism.
   * @param delta the increment 
   */
  virtual void Extend (size_t ) {};

  /**
   * copy an entry. 
   * This method has to be overwritten in child-classes to ensure a proper 
   * working of the cleaning mechanism.
   * @param src the index of the source entry
   * @param dest the index of the destination entry 
   */
  virtual void CopyEntry (size_t , size_t ) {};

  /**
   * notify a change of the list structure.
   * This method has to be called, whenever the list structure has been altered.
   * It results in an increased state counter.
   */
  void NotifyChange() { ++state_counter; };

  /// lock the mutex
  //  void LockMutex() { pthread_mutex_lock(&(master->mutex)); };

  /// release the mutex
  //  void UnlockMutex() { pthread_mutex_unlock(&(master->mutex)); };

public:
  /**
   * get the master of this List.
   * @return a pointer to the master List
   */
  List* Master() { return master; };

  /**
   * check if this is a master List.
   * @return true if it is a master List
   */
  bool IsMaster() { return this == master; };

  /**
   * get the maximal number of entries.
   * @ return the maximal number of entries
   */
  size_t MaxNumEntries() const { return max_num_entries; };

  /**
   * get the current number of entries.
   * @return the current number of entries
   */
  size_t NumEntries() const;

  /**
   * get the dynamical increment.
   * @return the dynamical increment
   */
  size_t DeltaEntries() const { return delta_entries; };

  /**
   * set the dynamical increment.
   * @param new_delta the new increment
   */
  void SetDelta(size_t new_delta);

  /**
   * constructor.
   * Creates an uninitialized List.
   */ 
  List ();

  /**
   * constructor. 
   * Creates a new List as master and initializes it.
   * @param mne initial maximum number of entries
   * @param delta initial increment of entries 
   */
  List (size_t mne, size_t delta);

  /**
   * constructor. 
   * Creates a new list which is linked to an existing master and
   * initializes it.
   * @param a_master the master List to link the new List to
   */
  List (List* a_master);

  /// copy constructor.
  List(const List &other);

  /// destructor.
  virtual ~List ();

  /**
   * get the active flag for an entry.
   * @param i index of the entry
   * @return the active flag
   */
  bool IsActive(size_t i) const {
    return master->active[i]; 
  };

  /**
   * get the number of active entries. 
   * @return the number of active entries
   */
  size_t NumActiveEntries () const;

  /**
   * add an entry at the end of the list.
   * @return the index of the new entry 
   */
  size_t AddEntry () { return Alloc(1); };

  /**
   * add an entry at the end of the list.
   * @return the index of the new entry 
   */
  void AddUntil (const size_t address) {
      while(NumEntries() <= address) {
	  AddEntry();
      };
  }

  /**
   * allocate a number of new entries.
   * @param n the number of new entries
   * @return the index of the first new entry
   */

  size_t Alloc(size_t n);

  /** 
   * find out which would be the first entry
   * if a block would be allocated.
   * @param n the number of new entries
   * @return the index of the first new entry
   */
  size_t TryAlloc(size_t n);// { return num_entries; };

  /**
   * delete an entry
   * @param i the index of the entry to be deleted 
   */
  void DelEntry (size_t i);

  /**
   * delete a number of continuous entries.
   * @param i the first entry to delete
   * @param n the number of entries to be deleted
   */
  void Delete(size_t i, size_t n) {
    while (n > 0) {
      DelEntry(i+n-1);
      n--;
    };
  };

  /**
   * delete the data of an entry.
   * This method has to be overwritten if a derived class
   * allocates memory for new entries.
   * @param i the index of the entry to be deleted 
   */
  virtual void DeleteData(size_t ) {};
    
  /// delete all entries.
  void DelAll ();

  /** 
   * find the next active entry. 
   * This method is a bit slow right now, but that
   * will be improved in the future. So please stick to this and do not try
   * to write a fast workaround.
   * @param i the index of the entry whose following entry has to be found
   * @return the index of the following entry 
   */
  size_t next_idx(size_t i) const;
  /**
   * find the previous active entry. 
   * This method is a bit slow right now, but that
   * will be improved in the future. So please stick to this and do not try
   * to write a fast workaround.
   * @param i the index of the entry whose following entry has to be found
   * @return the index of the following entry 
   */
  size_t prev_idx(size_t i) const;

  /** 
   * find the first active entry.
   * @return the index of the first active entry 
   */
  size_t begin_idx() const;

  /**
   * find the last active enty.
   * @return the index of the last active entry 
   */
  size_t last_idx() const;

  /**
   * specifies an index after the last entry.
   * @return an index after the last one
   */
  //size_t& end_idx() const; und size_t end_idx() const;
  //verhalten sich gleich size_t const end_idx() const ebenso
  // und der auch const size_t end_idx() const;
  size_t end_idx() const;
  

  /// create the offset-list for the cleaning process.
  void CreateOffsetList ();

  /// delete the offset list. 
  void DeleteOffsetList () { 
    delete [] offset; 
    offset = NULL;
  };

  /**
   * get the offset of an entry.
   * @param i the index of the entry
   * @return the offset of the entry 
   */
  size_t Offset (size_t i);

  /**
   * find out if the list has deleted entries or if it is clean.
   * @return true if it is clean
   */ 
  bool IsClean() const { return master->is_clean; };

  /**
   * clean the list.
   * This method eliminates all holes from the list and resets num_entries.
   * After CleanUp num_entries will be identical with the number of active entries. 
   */
  virtual void CleanUp ();

  /**
   * links this List to another List.
   * Please note that this can only be done for a 'standalone' List.
   * As soon as the List is linked to another list this method will
   * abort with an error message.
   * @param linked_list the List to link to
   */
  void Link (List *linked_list);

  /**
   * Check if a List is linked to annother one
   * @param other_list the list for which a link is checked
   * @param is_master_of_other in case of a link, indicate wether this is
   *        master of other_list (return reference).
   * @return bool indicating if "this" is anyhow linked to other_list
   */
    bool IsLinkedTo(List *other_list,
		    bool& is_master_of_other);

  /**
   * copy an entry. 
   * @param src the index of the source entry
   * @param dest the index of the destination entry 
   */
  virtual void Copy (size_t src, size_t dest);

  /**
   * check if an index is valid.
   * This method check if an index is valid in this List.
   * If the index is invalid an InvalidIndex_error will be
   * thrown.
   * @param i the index to check
   */
  void CheckIndex(size_t i) { 
    if (!IsActive(i)) {
      cerr << "index " << i << " is invalid" << endl;
      throw InvalidIndex_error(i); 
    };
  };

  /**
   * check if another List has the same structure than this one.
   * @param other_list the List to compare
   * @return true if the structure is the same
   */
  bool HasSameStructure(List *other_list);

  /**
   * Get the state counter of the list.
   * The state counter can be used to determine if the structure has changed
   * (e.g. entries have been allocated or deleted).
   * Whenever the List changes the state counter will be increased.
   * @return the state counter
   */
  state_t State() { return master->state_counter; };

  void operator=(const List &other);

    /**
     * swap two entries.
     * This method has to be overwritten in child-classes if Swap-processes will be applied
     * The effect is, that data at address i1 and i2 are exchanged
     * @param i1 the one index
     * @param i2 the other index
     */
    virtual void Swap(size_t i1, size_t i2){};

    /** swap two entries for all linked List
     *  @param i1 the one index
     *  @param i2 the other index
     */ 
    void SwapAll(size_t& i1, size_t& i2) {
	if (IsMaster()) {
	    Swap(i1, i2);
            /// @todo why unsigned int and not size_t ?
            for (unsigned int i = 0; i < num_clients; i++) {
		client[i]->Swap(i1, i2);
	    };
	} else {
	    Master()->SwapAll(i1, i2);
	};
    };

// ATTENTION: virtual methods are compiled for ALL TEMPLATE ARGUMENTS
// This causes problems due to the > operator. Use non virtual method instead
//    /**
//       Sorting mechanism for Lists. Virtual entry point only.
//    */
//    virtual void SortList() {
//	cout << "Error: Attempt to sort a List (method unknown)" << endl;
//	exit(EXIT_FAILURE);
//    };
};

inline size_t List::end_idx() const
{ 
  return master->after_last_entry; 
}

inline size_t List::begin_idx() const
{
  size_t i = 0;
  while ((i < end_idx()) && !IsActive(i)) ++i;
  return i;
}

inline size_t List::next_idx(size_t i) const 
{
  do { i++; } while ((i < end_idx()) && !IsActive(i));
  return i;

//     if(IsClean()) {
// 	i++;
// 	return i;
//     } else {
// 	do { i++; } while ((i < end_idx()) && !IsActive(i));
// 	return i;
//     };
}

inline size_t List::prev_idx(size_t i) const 
{
  if (i == 0) throw NotFound_error("List::prev");
  do { i--; } while ((i > 0) && !IsActive(i));
  if (!IsActive(i)) {
    throw NotFound_error("List::prev");
  };
  return i;
}

inline size_t List::last_idx() const 
{ 
  if (end_idx() == 0) {
    throw NotFound_error("List::last");
  };
  return end_idx() - 1; 
}

inline size_t List::NumEntries() const
{
  return end_idx();
}

/** 
 * a loop over all active entries of a List.
 * This macro expands to for(size_t I = L begin(); I < L end(); I = L next(I)).
 * @param I the index variable
 * @param L the List, including the member-operator (for example "this->") 
 */
#define FORALL(I,L) for (size_t I = L begin_idx(); I < L end_idx(); I = L next_idx(I))

/**
 * if IsClean() was done, the List starts at 0 and there are no inactive entries between 
 * 0 and end_idx()
 * a loop over all active entries of a clean List.
 * if the List is not clean, than the program exits!
 * @param I the index variable
 * @param L the List, including the member-operator (for example "this->") 
 */
#define FORALLCLEAN(I,L) if(!(L IsClean())){cout<<cout<<"List not clean"<<endl;exit(0);} else for (size_t I = 0; I < L end_idx(); I ++)


END_MOUSE

#endif












