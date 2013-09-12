#ifndef COMPRESSIBLESOLIDWALLBOBC_H
#define COMPRESSIBLESOLIDWALLBOBC_H

class CompressibleSolidWallBOBC;

#include "blockobjectbc.h"

class CompressibleSolidWallBOBC : public BlockObjectBC
{
public:
    CompressibleSolidWallBOBC(size_t field);

    CompressibleSolidWallBOBC(size_t field, BlockObject* block_object);

    virtual void operator()();

};

#endif // COMPRESSIBLESOLIDWALLBOBC_H
