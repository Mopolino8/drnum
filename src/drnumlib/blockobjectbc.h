#ifndef BLOCKOBJECTBC_H
#define BLOCKOBJECTBC_H

class BlockObjectBC;

#include "genericoperation.h"
#include "blockobject.h"

class BlockObjectBC : public GenericOperation
{

protected: // attributes
  BlockObject* m_BlockObject; // the block object to work on


public:

  BlockObjectBC();

  void setBlockObject(BlockObject* block_object) {m_BlockObject = block_object;}

  virtual void operator()() = 0;

};

#endif // BLOCKOBJECTBC_H
