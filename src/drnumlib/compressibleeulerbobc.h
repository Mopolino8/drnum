#ifndef COMPRESSIBLEEULERBOBC_H
#define COMPRESSIBLEEULERBOBC_H

class CompressibleEulerBOBC;

#include "blockobjectbc.h"

class CompressibleEulerBOBC : public BlockObjectBC
{
public:
    CompressibleEulerBOBC();
    virtual void operator()() {BUG;}
};

#endif // COMPRESSIBLEEULERBOBC_H
