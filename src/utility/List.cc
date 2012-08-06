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
// This program is distributed in the hope thaolst it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// !!

#include "List.hh"

#include <iostream>
#include <stdlib.h>

BEGIN_MOUSE

void List::InitList(size_t a_max_num_entries, size_t a_delta_entries)
{
    if (initialized) {
    cout << "fatal error: attempt to initialize a List twice" << endl;
    exit (EXIT_FAILURE);
  };
  // last_idx() bringt fehler
  max_num_entries = a_max_num_entries;
  delta_entries   = a_delta_entries;
  //ATTENTION: Changed, July/01/04
  if(max_num_entries < 1) {
      max_num_entries = 1;
  };
  if(delta_entries < 1) {
      delta_entries = 1;
  };
  //END
  offset = NULL;
  if (master == this) {

    active = new bool[max_num_entries];
    for (size_t i_entry = 0; i_entry < max_num_entries; i_entry++) {
      active[i_entry] = false;
    };
    block_size = new size_t[max_num_entries];
    UpdateBlockSize();
    after_last_entry = 0;
    first_block = 0;

    // create pthread-mutex
    //    pthread_mutex_init(&mutex, NULL);

  } else {
    active = NULL;
    block_size = NULL;
  };
  client = NULL;
  num_clients = 0;
  initialized = true;
}

// kann evtl geloescht werden
List::List()
{
  client = NULL;
  master = this;
  initialized = false;
  num_clients = 0;
  state_counter = 0;
  offset = NULL;
  active = NULL;
  block_size = NULL;
  after_last_entry = 0;
}

List::List(size_t a_max_num_entries, size_t a_delta_entries)
{
  is_clean = true;
  initialized = false;
  master = this;
  InitList (a_max_num_entries, a_delta_entries);
  num_clients = 0;
}

List::List(List *linked_list) 
{
  if(linked_list == NULL) {
      cout << " atempt to initialize List with NULL master" << endl;
      cout << " List::List(List *linked_list) " << endl;
      exit(EXIT_FAILURE);
  };
  client = NULL;
  initialized = false;
  master = linked_list->master;
  InitList (master->MaxNumEntries(), master->DeltaEntries());
  master->AddClient(this);
  num_clients = 0;
}

List::List(const List &other)
{
  if(&other == NULL) {
      cout << " atempt to initialize List with NULL master" << endl;
      cout << " List::List(const List &other) " << endl;
      exit(EXIT_FAILURE);
  };
  client = NULL;
  master = this;
  initialized = false;
  num_clients = 0;
  state_counter = 0;
  offset = NULL;
  active = NULL;
  block_size = NULL;
  after_last_entry = 0;
  operator=(other);
}

void List::operator=(const List &other)
{
  if(&other == NULL) {
      cout << " atempt to initialize List with NULL master" << endl;
      cout << " List::operator=(const List &other) " << endl;
      exit(EXIT_FAILURE);
  };
  if (initialized) {
    cerr << "fatal error: operator=() can only be udes for uninitialized Lists." << endl;
    exit(EXIT_FAILURE);
  };
  if (other.master == &other) master = this;
  else master = other.master;
  if (other.initialized) {
    InitList(other.max_num_entries, other.delta_entries);
    // Fake that we are extending the list.
    // This is done to make sure that necessary memory will be allocated in derived classes
    size_t mne = max_num_entries;
    max_num_entries = 0;
    Extend(mne);
    max_num_entries = mne;
  };
  if (this != master) {
    master->AddClient(this);
  } else if (initialized) {
    after_last_entry = other.after_last_entry;
    for (size_t i = 0; i < max_num_entries; ++i) active[i] = other.active[i];
    UpdateBlockSize();
  };
}

void List::Link (List *linked_list) 
{
  if(linked_list == NULL) {
      cout << " atempt to initialize List with NULL master" << endl;
      cout << " void List::Link (List *linked_list) " << endl;
      exit(EXIT_FAILURE);
  };
  if (this != master) {
    cout << "This 'List' is already linked to another 'List'!" << endl;
    exit(EXIT_FAILURE);
  };
  if (num_clients != 0) {
    cout << "This 'List' is already a master for someone!" << endl;
    exit(EXIT_FAILURE);
  };
  client = NULL;
  master = linked_list->master;
  DelAll();
  InitList (master->MaxNumEntries(), master->DeltaEntries());

  // Fake that we are extending the list.
  // This is done to make sure that necessary memory will be allocated in derived classes
  size_t mne = max_num_entries;
  max_num_entries = 0;
  Extend(mne);
  max_num_entries = mne;

  FORALL(i, this->) InitEntry(i);
  master->AddClient(this);
}

bool List::IsLinkedTo (List *other_list,
		       bool& is_master_of_other)
{
    if(this) {
	if(other_list == NULL) {
	    is_master_of_other = false;
	    return false;
	};

	// "this" and other_list are linked toigether, if the master is the same
	bool linked;
	linked = (other_list->master == master);

	// check, if "this is master of other_list
	is_master_of_other = (linked && (this == master));

	return linked;
    };
    return false;
}

void List::AddClient (List *client_to_add) 
{
  if (this != master) {
    cout << "fatal error trying to link to a non master list!" << endl;
    exit(EXIT_FAILURE);
  };
  List **new_client;
  new_client = new List* [num_clients + 1];
  for (unsigned int i = 0; i < num_clients; i++) new_client[i] = client[i];
  new_client[num_clients] = client_to_add;
  num_clients++;
  delete [] client;
  client = new_client;
  /*
  while (client_to_add->num_entries < num_entries) {
    client_to_add->InitEntry(client_to_add->num_entries);
    client_to_add->num_entries++;
  };
  */
}

void List::DelClient(List *client_to_del){
  if (this != master) {
    cout << "fatal error trying to delete a client from a non master list!" << endl;
    exit(EXIT_FAILURE);
  };
  List **new_client;
  int off = 0;
  for (unsigned int i = 0; i < num_clients - off; i++) {
    if (client[i] == client_to_del) off += 1;
    if((i + off) < num_clients) {   //!!! fixed invalid read access 10/28/02 R.
      client[i] = client[i + off];  //!!! minor error
    };
  };
  new_client = new List* [num_clients - off];
  for (unsigned int i = 0; i < num_clients- off; i++) new_client[i] = client[i];
  num_clients -= off;
  delete [] client;
  client = new_client;
  NotifyChange();
}

List::~List () 
{
    if (this != master) {
	master->DelClient(this);
	delete [] offset; // ATTENTION: not sure if needed, no master
	delete [] client; // ATTENTION: not sure if needed, no master
	delete [] active; // ATTENTION: not sure if needed, no master
    } else {
	if(num_clients >= 1) { // if clients exist
	    // make your first client be the new master
	    List* new_master = client[0];
	    new_master->master = new_master;
	    // take out previous client[0], now new_master, from client list, as it will be the new master
	    DelClient(new_master);
	    // tell all remaining clients, who's the new master
	    for (size_t i = 0; i < num_clients; i++) {
		client[i]->master = new_master;
	    };
	    // let the new master take over master-specific arrays and data
	    new_master->block_size = block_size;
	    new_master->offset = offset;
	    new_master->client = client;
	    new_master->active = active;
	    new_master->first_block = first_block;
	    new_master->after_last_entry = after_last_entry;
	    new_master->delta_entries = delta_entries;
	    new_master->max_num_entries = max_num_entries;
	    new_master->max_offset = max_offset;
	    new_master->num_clients = num_clients;
	    new_master->is_clean = is_clean;
	    new_master->block_size_updated = block_size_updated;
	    new_master->state_counter = state_counter;
	    new_master->initialized = initialized;
	} else { // no clients, delete all master-specific arrays
	    delete [] block_size;
	    delete [] offset;
	    delete [] client;
	    delete [] active;
	};
    };
}

size_t List::NumActiveEntries () const
{
  size_t num_actives = 0;
  for (size_t i_entry = begin_idx(); i_entry < end_idx(); i_entry = next_idx(i_entry)) {
    num_actives++;
  };
  return num_actives;
}

void List::ExtendList (size_t delta) 
{
  if (this != master) {
    cout << "fatal error in ExtendList this != master." << endl;
    exit(EXIT_FAILURE);
  };
  bool *new_active;
  Extend (delta);
  max_num_entries += delta;
  new_active = new bool[max_num_entries];
  for (size_t i = 0; i < max_num_entries; i++) {
    if (i < max_num_entries - delta) {
      new_active[i] = active[i];
    } else {
      new_active[i] = false;
    }
  };
  delete [] active;
  active = new_active;
  for (unsigned int j = 0; j < num_clients; j++) {
    client[j]->Extend(delta);
    client[j]->max_num_entries += delta;
  };
  delete [] block_size;
  block_size = new size_t[max_num_entries];
  UpdateBlockSize();
}

void List::DelEntry (size_t i_entry)
{
  if (!IsActive(i_entry)) throw InvalidIndex_error(i_entry);
  for (unsigned int i_client = 0; i_client < num_clients; i_client++) {
    client[i_client]->DeleteData(i_entry);
  };
  DeleteData(i_entry);
  master->active[i_entry] = false;
  master->block_size[i_entry] = 1;
  is_clean = false;
  NotifyChange();
}

void List::DelAll () {
  if (initialized) {
    for (size_t i_entry = 0; i_entry < max_num_entries; i_entry++) {
      if (IsActive(i_entry)) DelEntry (i_entry);
    };
    CleanUp();
    NotifyChange();
  };
}

void List::CreateOffsetList () {
  if (master == this) {
    max_offset = 0;
    DeleteOffsetList();
    offset = new size_t [end_idx()];
    for (size_t i_entry = 0; i_entry < end_idx(); i_entry++) {
      offset[i_entry] = max_offset;
      if (!IsActive(i_entry)) max_offset++;
    };
  } else {
    master->CreateOffsetList();
  };
}

size_t List::Offset (size_t i_entry) {
  if (master->offset == NULL) {
    cerr << "error: no offset information available" << endl;
    cerr << "       at this point in the program\n" << endl;
    exit (EXIT_FAILURE);
  } else {
    return master->offset[i_entry];
  };
  return 0; // dummy
}

void List::Copy(size_t src, size_t dest)
{
  if (IsMaster()) {
    CopyEntry(src, dest);
    for (unsigned int i = 0; i < num_clients; i++) {
      client[i]->CopyEntry(src, dest);
    };
  } else {
    Master()->Copy(src, dest);
  };
}

void List::CleanUp () {
  size_t i_entry;
  if (master == this) {
    CreateOffsetList();
    for (i_entry = 0; i_entry < end_idx(); i_entry++) {
      if (Offset(i_entry) != 0) {
	Copy(i_entry, i_entry - Offset(i_entry));
      };
      if (i_entry < end_idx() - max_offset) {
	active[i_entry] = true;
      } else {
	active[i_entry] = false;
      };
    };
    UpdateBlockSize();
    try {
      //after_last_entry = prev_idx(max_num_entries) + 1;
      size_t i = max_num_entries;
      if (i == 0) {
	throw NotFound_error("List::prev");
      };
      do {
	i--;
      } while ((i > 0) && !IsActive(i));
      if (!IsActive(i)) {
	throw NotFound_error("List::prev");
      };
      after_last_entry = i + 1;
    } catch (NotFound_error) {
      //DBG_CP;
      after_last_entry = 0;
      //DBG_CP;
    };
    is_clean = true;
  } else {
    master->CleanUp();
  };
}

void List::SetDelta(size_t new_delta)
{
  if (this != master) {
    cerr << "cannot modify the delta value for client lists" << endl;
    exit(EXIT_FAILURE);
  };
  delta_entries = new_delta;
}

void List::UpdateBlockSize()
{
  if (max_num_entries != 0) {
    block_size[max_num_entries - 1] = 1;
    for (size_t i = max_num_entries - 1; i > 0; i--) {
      if ((active[i] && !active[i-1]) || (!active[i] && active[i-1])) {
	block_size[i-1] = 1;
      } else {
	block_size[i-1] = block_size[i] + 1;
      };
    };
  };
  block_size_updated = true;
  first_block = 0;
}

size_t List::TryAlloc(size_t n)
{
  if (IsMaster()) 
    {
      size_t new_entry = first_block;
      bool found = false;
      do 
	{
	  if (!active[new_entry] && (block_size[new_entry] >= n)) 
	    {
	      found = true;
	    } 
	  else 
	    {
	      new_entry += block_size[new_entry];
	    };
	} while ((new_entry < max_num_entries) && !found);
      if (found) 
	{
	  return new_entry;
	} 
      else 
	{
	  if (block_size_updated) 
	    {
	      ExtendList(delta_entries);
	      block_size_updated = false;
	    } 
	  else 
	    {
	      UpdateBlockSize();
	    };
	  return TryAlloc(n);
	};
      //never reached : changed by kaiser 12.01.06
      //return after_last_entry;
    } 
  else 
    {
      return master->TryAlloc(n);
    };
}

size_t List::Alloc(size_t n)
{
  if (IsMaster()) {
    if (!initialized) {
      size_t d = n/10;
      if (d == 0) d = 1;
      InitList(0,d);
      ExtendList(n);
    };
    size_t new_entry = TryAlloc(n);
    for (size_t i = new_entry; i < new_entry + n; i++) {
      block_size[i] = n - i + new_entry;
      active[i] = true;
      InitEntry(i);
      for (unsigned int j = 0; j < num_clients; j++) {
	client[j]->InitEntry(i);
      };
    };
    if (new_entry+1 > after_last_entry) after_last_entry = new_entry + n;
    first_block = new_entry;
    NotifyChange();
    return new_entry;
  } else {
    return master->Alloc(n);
  };
}

bool List::HasSameStructure(List *other_list)
{
  if (begin_idx() != other_list->begin_idx()) return false;
  if (end_idx() != other_list->end_idx()) return false;
  if (MaxNumEntries() != other_list->MaxNumEntries()) return false;
  for (size_t i = 0; i < MaxNumEntries(); i++) {
    if (IsActive(i) != other_list->IsActive(i)) return false;
  };
  return true;
}


END_MOUSE





