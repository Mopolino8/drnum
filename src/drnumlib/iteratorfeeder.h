#ifndef ITERATORFEEDER_H
#define ITERATORFEEDER_H

//#include <cstddef>
//#include <string.h>
//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>

#include "blockcfd.h"
#include "iterators/patchiterator.h"
#include "patchgrid.h"

#include <string>
using namespace std;

//class IteratorFeeder;

class IteratorFeeder
{

  vector<PatchIterator*> m_Iterators;

protected: // methods

public: // methods

  /**
   * Add an Iterator
   * @param
   */
  void addIterator(PatchIterator* iterator);

  void feed(PatchGrid& patch_grid);

};

#endif // ITERATORFEEDER_H
