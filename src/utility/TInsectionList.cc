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

template <class T>
TInsectionList<T>::TInsectionList (List *a_main_list,
				   size_t a_default_size, size_t a_increment)
    : List(a_main_list)
{
    // Keep pointer to main_list
    main_list = a_main_list;

    // size parameters for item_list to build
    default_size = a_default_size;
    increment = a_increment;
    if(default_size == 0) {
	default_size = main_list->NumEntries();
    };
    if(increment == 0) {
	increment = main_list->NumEntries();
    };

    // Build the start list
    start = new TList<size_t>(a_main_list);

    // Set item_list and item_main_add to not available
    item_list = NULL;
    item_main_add = NULL;
}

template <class T>
void TInsectionList<T>::Initialize()
{
    // Destroy previously allocated data, if existent
    //.. item_list
    if(item_list) {
	delete item_list;
    };
    if(item_main_add) {  // should never be true, since item_main_add is linked to item_list
	delete item_main_add;
    };
    item_list = new TList<T>(default_size, increment);
    item_main_add = new TList<size_t>(item_list);
}

template <class T>
void TInsectionList<T>::Finalize(bool sparse)
{
    // Resort the item_list and along with it also item_main_add. Sort by addresses of item_main_add.
    // Do the following:
    // 1) Count number of occurencies for the main_list_address. Use "start" as memory.
    // 2) Compute start addresses.
    // 3) Go through list and resort (new according to Rainer):
    //    * Use while-loop to go totally once through the item_list, index old_add.
    //    * Find new_add for the given item_list-entry.
    //     * Do not swap, but increment old_add, if item is in correct range.
    //     * Else, swap new_add and old_add contents.
    // 4) If sparse == true, item_main_add is deleted after list is sorted.
    // 5) Compute maximum number of items for any of the main_list entries.
    
    // 1) Count number of occurencies for the main_list_address.
    //    REMARK: preliminary use of start as counter
    FORALL(i_m, main_list->) {
	start->At(i_m) = 0;
    };
    FORALL(i, item_main_add->){
	start->At(item_main_add->At(i))++;  // note useage as counter!!!
    };

    // 2) Compute start address:
    //    REMARK: meaning of "start" changes from counter to actual start address
    size_t add_count = 0;
    FORALL(i_m, main_list->){
	size_t add_count_next = add_count + start->At(i_m);
	start->At(i_m) = add_count;
	add_count = add_count_next;
    };

    // 3) Go through list and resort (new according to Rainer):
    // Start sorting at the beginning of the list
    size_t old_add = 0;
    // Use while-loop to go totally once through the item_list.
    while(old_add < item_list->NumEntries()) {
	bool do_swap = true;
	size_t main_list_add = item_main_add->At(old_add);
	size_t new_add = start->At(main_list_add);
	// Run over all items, that have the same main_list address
	while(item_main_add->At(new_add) == main_list_add) {
	    if(new_add == old_add) {
		do_swap = false;
		break;
	    };
	    new_add++;
	};
//not needed if(new_add == old_add) {
//not needed	    do_swap = false;
//not needed	};
	if(do_swap) {
	    // exchange list entries of new_add and old_add
	    item_list->SwapAll(new_add, old_add);
	} else {
	    // nothing else to do with old_add, go on further
	    old_add++;
	};
    };

// Debug only
//    WriteActual();

    // Check, if sorting was successfull.
    // ATTENTION: Discard this test after test-phase of this class
    size_t start_address;
    size_t end_address;
    FORALL(i_m, main_list->) {
	start_address = start->At(i_m);
	if(i_m < (main_list->NumEntries()-1)) {
	    end_address = start->At(i_m + 1);
	} else {
	    end_address = item_list->NumEntries();
	};
	for(size_t i=start_address; i<end_address; i++) {
	    if(item_main_add->At(i) != i_m) {
		cout << "ERROR, TInsectionList, Wrongly sorted item_list" << endl;
		exit(EXIT_FAILURE);
	    };
	};
    };

    // 4) If a sparse == true, item_main_add is deleted after list is sorted.
    if(sparse) {
	delete item_main_add;
    };

    // 5) Compute maximum number of items in second dimension.
    ComputeMaxNumItemsSecondDim();
}

template <class T>
void TInsectionList<T>::ComputeMaxNumItemsSecondDim()
{
    max_num_items_second_dim = 0;
    FORALL(i_m, main_list->) {
      size_t num_items_second_dim;
      if(i_m < (main_list->NumEntries()-1)) {
	num_items_second_dim = start->At(i_m + 1) - start->At(i_m);
	//	end_address = start->At(i_m + 1);
      } else {
	num_items_second_dim = item_list->NumEntries() - start->At(i_m);
      };
      if(num_items_second_dim > max_num_items_second_dim) {
	max_num_items_second_dim = num_items_second_dim;
      }
    }
}

template <class T>
void TInsectionList<T>::MakeUnique()
{
    // Eliminate all multiple items in second dimension
    // Do the following:
    // Go through main_list (first dimension, index i_m), go through item_list in range of i_m
    // and eliminate duplicates by copying uniques back to consecutive positions.
    // Array start ist adjusted by incrementing down to new starting addresses of ranges.
    // Delete obsolete rest of array at the end to free memory space (tight fit).
    //
    // ATTENTION: O(N2nd^2): Quadratic bahaviour in second dimension.
    //            Use a sort algorithm and then an "=" based elimination algorithm later.
    //            This avoids quadratic behaviour in second dimension.
    //
    // Prime
    //.. Set a running address to copy unique entries
    size_t j_target = 0;
    //.. Loop for 1st dimension
    FORALL(i_m, main_list->){
	//.. Prime
	//.... shifted new start address for items in range of i_m
	size_t start_newadd = j_target;
	//.... range of addresses to check, before shifting down
	size_t start_oldadd = start->At(i_m);
	size_t after_end_oldadd;
	if(i_m < NumEntries()-1) {
	    after_end_oldadd = start->At(i_m+1);
	} else {
	    after_end_oldadd = item_list->NumEntries();
	};
	//.. Loop for items to test (starts at first item, to ensure copying to
	//   shifted new location, even though first item in range never has a previous same)
	for(size_t j = start_oldadd; j < after_end_oldadd; j++) {
	    //.... Check, if item_list->At(j) is equal to any of the previous known
	    bool duplicate = false;
	    for(size_t j_known = start_newadd; j_known < j_target; j_known++) {
		if(item_list->At(j) == item_list->At(j_known)) {
		    duplicate = true;
		    break;
		};
	    };
	    if(!duplicate) { // item not occured yet, copy onto new position
		//...... Copy entry to new positions. Applies to all linked lists
		item_list->Copy(j, j_target); // will copy linked item_main_list too
		j_target++;
	    };
	};
	//.. Adjust start->At(i_m)
	start->At(i_m) = start_newadd;
    };
    //
    // 2) Delete rest of array item_list at the end, that is obsolete after shifting items back
    //    to unique sets
    item_list->Delete(j_target, (item_list->NumEntries() - j_target));
    item_list->CleanUp();

    // 3) Compute maximum number of items in second dimension.
    ComputeMaxNumItemsSecondDim();
}

// template <class T>
// void TInsectionList<T>::MakeUniqueSorted()
// {
//     // Sort and eliminate all multiple items in second dimension
//     // Do the following:
//     // Go through main_list (first dimension, index i_m), sort item_list in range of i_m,
//     // go through sorted item_list in range of i_m and eliminate duplicates by copying uniques back to
//     // consecutive positions.
//     // Array start ist adjusted by incrementing down to new starting addresses of ranges.
//     // Delete obsolete rest of array at the end to free memory space (tight fit).
//     //
//     // ATTENTION: O(N2nd * ln(N2nd)): Logarithmic bahaviour (following the sort implementation from
//     //            STL-lirary) in second dimension.
//     //
//     // Prime
//     //.. Set a running address to copy unique entries
//     size_t j_target = 0;
//     //.. Loop for 1st dimension
//     FORALL(i_m, main_list->){
// 	//.. Prime
// 	//.... shifted new start address for items in range of i_m
// 	size_t start_newadd = j_target;
// 	//.... range of addresses to check, before shifting down
// 	size_t start_oldadd = start->At(i_m);
// 	size_t after_end_oldadd;
// 	if(i_m < NumEntries()) {
// 	    after_end_oldadd = start-At(i_m+1);
// 	} else {
// 	    after_end_oldadd = item_list->NumEntries();
// 	};
// 	//.. Sort items in range of i_m
// 	//.... Transfer to temporary vector
// 	vector<pair<T, size_t> > v_help;
// 	v_help.clear();
// 	for(size_t j = start_oldadd; j < after_end_oldadd; j++) {
// 	    pair<T, size_t> pp;
// 	    pp.first = item_list->At(j);
// 	    pp.second = j;
// 	    v_help.push_back(pp);
// 	};
// 	//.... Sort v_help
// 	v_help.sort(v_help.begin(), v_help.end());
// 	//.... Get through sorted v_help and eliminate non uniques. This operation is of lin. complexity
// 	//     ATTENTION: Cant use pre-defined STL unique, since it reacts on pair::second as well. To
// 	//                use STL unique, must write a comparison op. Direct code written here.
// 	bool first = true;
// 	for(vector<size_t>::iterator k_v_help = v_help.begin(); k_v_help != v_help.end(); k_v_help++) {

// // HIER GEHTS WEITER:
// // Problem: gelinkte Listen zu item_list, vor allem item_main_add werden müssen ggf. umkopiert werden.
// //          Dabei darf keine Überschreibung kleinerer Adressen erfolgen. Es muss mit SwapAll statt
// //          Copy gearbeitet werden. Erst alle Swappen, dann löschen.

// 	    T item_j = *v_help.first;
// 	    size_t j = *v_help.second;
// 	    if(first){   // do not ask for uniquenes of first item in range
// 		item_list->Copy(j, j_target);
// 		j_target++;
// 		first = false;
// 	    };
// 	    if(*v_help.first != *(v_help-1).first) {  // item is not equal to previous
		


// 	vector<T>::iterator new_end = v_help.unique(v_help.begin(), v_help.end());
// 	//.... Transfer v_help in unque range back to item_list
// 	for(vector<size_t>::iterator k_v_help = v_help.begin(); k_v_help != v_help.end(); k_v_help++) {
	    


// 	//.. Loop for items to test (starts at first element, to ensure copying to
// 	//   shifted new location, even though first item in range never has a previous same)
// 	for(size_t j = start_oldadd; j < after_end_oldadd; j++) {
// 	    //.... Check, if item_list->At(j) is equal to any of the previous known
// 	    bool duplicate = false;
// 	    for(size_t j_known = start_newadd; j_known < j_target; j_known++) {
// 		if(item_list->At(j) == item_list->At(j_known)) {
// 		    duplicate = true;
// 		    break;
// 		};
// 	    };
// 	    if(!duplicate) { // item not occured yet, copy onto new position
// 		//...... Copy entry to new positions. Applies to all linked lists
// 		item_list->Copy(j, j_target); // will copy linked item_main_list too
// 		j_target++;
// 	    };
// 	};
// 	//.. Adjust start->At(i_m)
// 	start->At(i_m) = start_newadd;
//     };
//     //
//     // 2) Delete rest of array item_list at the end, that is obsolete after shifting items back
//     //    to unique sets
//     item_list->Delete(j_target, (item_list->NumEntries() - j_target));
//     item_list->CleanUp();
//
// // 3) Compute maximum number of items in second dimension.
// ComputeMaxNumItemsSecondDim();
// }

template <class T>
void TInsectionList<T>::WriteActual()
{
    FORALL(i_m, main_list->){
	cout << start->At(i_m) << " ";
    };
    cout << endl;
    cout << endl;

//    FORALL(i, item_main_add->){
//	cout << item_list->At(i) << " ";
//    };
    cout << endl;
    FORALL(i, item_main_add->){
	cout << item_main_add->At(i) << " ";
    };
    cout << endl;
    cout << endl;
}


//     // 3a) Prime 1st unsorted element -> set on new position and help-store the element previously there.
//     TList<bool>* done;
//     done = new TList<bool>(item);
//     size_t old_add = 0;
//     T item_old_help = item->At(old_add);
//     T item_new_help;

//     while(total_flip_count < item->NumEntries()) {
// 	bool continue = true;
// 	while(continue) {
// 	    size_t new_add = new_from_old->At(old_add);
// 	    if(done->At(new_add)) {
// 		continue = false;
// 	    } else {
// 		item_new_help = item->At(new_add);
// 		item->At(new_add) = item_old_help;
// 		item_old_help = item_new_help;
// 		old_add = new_add;
// 		done->At(new_add) = true;
// 	    };
// 	};
//     };



//     // Neu nach Rainer:
//     // Start sorting at the beginning of the list
//     size_t old_add = 0;
//     // Check, if done
//     while(old_add < item->NumEntries()) {
// 	bool continue = true;
// 	main_list_add = item_main_add->At(old_add);
// 	new_add = start->At(main_list_add);
// 	// Run over all items, that have the same main_list address
// 	while(item_main_add->At(new_add) == main_list_add) {
// 	    if(new_add == old_add) {
// 		continue = false;
// 		new_add++;
// 	    };
// 	};
// 	if(new_add == old_add) {
// 	    continue = false;
// 	};
// 	if(continue) {
// 	    // exchange list entries of new_add and old_add
// 	    item->SwapAll(new_add, old_add);
// 	} else {
// 	    // nothing else to do with old_add, go on further
// 	    old_add++;
// 	};
//     };



// template <class T>
// void TInsectionList::Finalize()
// {
//     // get starting addresses in item_list
//     size_t current_address = 0;
//     FORALL(i, intermediate_count->) {
// 	start->At(i) = current_address;
// 	current_address+=intermediate_count->At(i);
//     };
//     size_t total_item_list_size = current_address;

//     // Total size of item_list is in current_address. Build item_list in exact size.
//     item_list = new TList<pair<size_t, T> >(total_item_list_size, 1);
//     item_list->Alloc(total_item_list_size);

//     // reset the counter
//     FORALL(i_c, intermediate_count->) {
// 	intermediate_count->At(i_c) = 0;
//     };

//     // loop through twin entries and set data in itemlist
//     for(vector<pair<size_t, T> >::iterator ite = intermediate_twin_entries->begin();
// 	ite != intermediate_twin_entries->end();
// 	ite++) {
// 	size_t main_list_address = ite->first;
// 	T item = ite->second;
// 	item_list->At(start->At(main_list_address)+intermediate_count->At(main_list_address))
// 	    = item;
// 	intermediate_count->At(main_list_address)++;
//     };

//     // delete obsolete data (hope he will actually do so)
//     delete intermediate_twin_entries;
//     delete intermediate_count;
// }

template <class T>
void TInsectionList<T>::Reset()
{
    // all reset functions also done in Initialize()
   TInsectionList<T>::Initalize();
}

template <class T>
TInsectionList<T>::~TInsectionList()
{
    if(start) {
	delete start;
    };
    if(item_list) {
	delete item_list; // this deletes item_main_add too
    };
}
