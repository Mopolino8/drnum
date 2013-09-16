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
#ifndef BLOCKOBJECTBC_H
#define BLOCKOBJECTBC_H

class BlockObjectBC;

#include "genericoperation.h"
#include "blockobject.h"

class BlockObjectBC : public GenericOperation
{

protected: // attributes
  BlockObject* m_BlockObject; // the block object to work on
  size_t m_Field;             // the variable field to work on


public:

  BlockObjectBC (size_t field);

  BlockObjectBC (size_t field, BlockObject* block_object);

  void setBlockObject (BlockObject* block_object) {m_BlockObject = block_object;}

  virtual void operator ()() = 0;

  void standardBlackCells0 ();
};

#endif // BLOCKOBJECTBC_H
