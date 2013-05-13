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

#ifndef TInsectionList_HH
#define TInsectionList_HH

#include <utility/namespace_mouse.hh>

BEGIN_MOUSE

template <class T> class TInsectionList;

END_MOUSE

#include <mouse_types.hh>
#include <TList.hh>
#include <vector>

BEGIN_MOUSE

/** A TInsectionList couples two lists (main_list and item_list) to store
 *  arbitrary ammounts of items in item_list for a single entry in main_list.
 *  Both lists are continuous in memory.
 *
 *  Data organisation:
 *    main_list:                         The list to which this TInsectionList is connected
 *    start: (linked to main_list):      The starting address in item_list
 *    item_list (non linked):            Sequence of items, type T, divided into branches
 *                                       for each entry of main_list.
 *    item_main_add (li. to item_list):  List linked to item_list, additionally containing
 *                                       the main_list address, the item belongs to.
 *
 *
 *  <pre>
 *  [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-] (item_list)
 *  [-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-][-] (item_main_add)
 *   ^        ^           ^        ^        ^     ^     ^        ^
 *   |   _____|           |        |        |     |     |        |
 *   |  |   ______________|        |        |     |     |        |
 *   |  |  |   ____________________|        |     |     |        |
 *   |  |  |  |   __________________________|     |     |        |
 *   |  |  |  |  |   _____________________________|     |        |
 *   |  |  |  |  |  |   ________________________________|        |
 *   |  |  |  |  |  |  |   ______________________________________|
 *   |  |  |  |  |  |  |  |
 *  [-][-][-][-][-][-][-][-] (start)
 *  </pre>

 *
 *  Automatic creation:
 *    ConsiderStore(.., ..) collects data pairs of type (T, main_list_address).
 *    Finalize() produce data lists as required and frees temporary memory.
 *    MakeUnique() produce data lists as Finalize, but eliminate multiple entries in 2nd dimension,
 *                 meaning all entries of item_list in range of start(nn) .. start(nn*1)-1 will
 *                 be made unique.
 */

template <class T>
class TInsectionList : public List
{
  protected:
  /** Pointer to a main_list list (no own data) */
  List* main_list;

  /** list of starting addresses in item_list (linked to main_list, persistent) */
  TList<size_t>* start;

  /** List to hold item data */
  TList<T>* item_list;

  /** List to hold main_list_address for items (linked to item_list). */
  TList<size_t>* item_main_add;

  size_t default_size;
  size_t increment;

  /** Maximum number of items for any of the main_list entries. This is the
     *  maximum size in the "second dimension". */
  size_t max_num_items_second_dim;

  /** Debug method to write start and item_main_add data.
     */
  void WriteActual();

  public:
  /** constructor
     *  @param a_main_list the main_list list to link to
     *  @param default_size the default size of the item_list to build
     *  @param increment the default incrementof the item_list to build
     *  NOTE: if default_size or increment are set to zero, the size of the
     *        main_list will be chosen.
     */
  TInsectionList(List *a_main_list, size_t default_size=0, size_t increment=0);

  /** Initialize data collecton. Allocates data lists needed.
     */
  void Initialize();

  /** Method to collect pairs of data.
     *  @param main_list_address the main_list address for which the item will be stored
     *  @param item_address the address to store
     */
  void ConsiderStore(const size_t& main_list_address, T& item) {
    // insert the pair in (non-sorted) item_list
    size_t item_address = item_list->AddEntry();
    item_list->At(item_address) = item;
    item_main_add->At(item_address) = main_list_address;
  }

  /** Method to Create data fields after all data pairs are stored.
     *  @param sparse true directs Finalize() to delete the item_main_add after
     *                sorting the item_list
     */
  void Finalize(bool sparse = false);

  /** Compute max number of items for any of the main_list entries. This is the
     *  maximum size in the "second dimension".
     */
  void ComputeMaxNumItemsSecondDim();

  /** Eliminate multiple entries in 2nd dimension: All entries of item_list in range of
     *  start(nn) .. start(nn*1)-1 will be made unique for each index nn in first dimension.
     */
  void MakeUnique();

  /** Sort and eliminate multiple entries in 2nd dimension: All entries of item_list in
     *  range of start(nn) .. start(nn*1)-1 will be sorted and made unique for each index nn
     *  in first dimension.
     */
  //    void MakeUniqueSorted();

  /** Method to Reset all stored data. After this, new allocations can be made.
     */
  void Reset();

  /** Get the number of entries in the second dimension.
     *  @param i the index in the main_list (first dimension)
     *  @return the number of items stored for i
     */
  size_t NumItems(size_t i) {
    if(i<(NumEntries()-1)) {
      return (start->At(i+1) - start->At(i));
    } else {
      return (item_list->NumEntries() - start->At(i));
    };
  }

  /** Get the total number of entires in the item_list.
     *  @return number of entires in the item_list
     */
  size_t NumItems() {
    return (item_list->NumEntries());
  }

  /** Access max number of items for any of the main_list entries. This is the
     *  maximum size in the "second dimension".
     *  @return number of entires in the item_list
     */
  size_t MaxNumItemsSecondDim() {
    return max_num_items_second_dim;
  }

  /** Data access method (read/write)
     *  @param i the address in the main_list (first dimension)
     *  @param k the index in the second dimension
     *  @return a reference to the entry
     */
  T& At(size_t& i, size_t& k) {
    return item_list->At(start->At(i)+k);
  }

  /** Get the address of an entry in the item_list
     *  @param i the address in the main_list (first address)
     *  @param k the index in the second dimension
     *  @return the address of the entry
     */
  size_t ItemIndex(size_t& i, size_t& k) {
    return start->At(i)+k;
  }

  /** Get a link point on the start array.
     *  @return pointer to start.
     */
  TList<size_t>* StartLinkPoint() {
    return start;
  }

  /** Get a link point on the item_list.
     *  @return pointer to item_list.
     */
  TList<T>* ItemLinkPoint() {
    return item_list;
  }

  /** Get a link point on the item_main_add.
     *  @return pointer to item_main_add.
     */
  TList<size_t>* ItemMainAddLinkPoint() {
    return item_main_add;
  }

  /** destructor */
  virtual ~TInsectionList ();
};

#include "TInsectionList.cc"
END_MOUSE

#endif
