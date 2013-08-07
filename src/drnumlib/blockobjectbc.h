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

  void standardBlackCells ();
};

#endif // BLOCKOBJECTBC_H
