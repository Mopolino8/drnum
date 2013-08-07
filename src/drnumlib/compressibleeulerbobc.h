#ifndef COMPRESSIBLEEULERBOBC_H
#define COMPRESSIBLEEULERBOBC_H

class CompressibleEulerBOBC;

#include "blockobjectbc.h"

class CompressibleEulerBOBC : public BlockObjectBC
{
public:
    CompressibleEulerBOBC(size_t field);

    CompressibleEulerBOBC(size_t field, BlockObject* block_object);

    virtual void operator()();
};

#endif // COMPRESSIBLEEULERBOBC_H
