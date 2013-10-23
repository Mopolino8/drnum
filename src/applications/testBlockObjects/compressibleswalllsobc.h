#ifndef COMPRESSIBLESWALLLSOBC_H
#define COMPRESSIBLESWALLLSOBC_H

class CompressibleSWallLSOBC;

#include "levelsetobjectbc.h"

class CompressibleSWallLSOBC : public LevelSetObjectBC
{
public:

    CompressibleSWallLSOBC(size_t field,
                           LevelSetObject* levelset_object,
                           size_t abuse_field);

    virtual void operator()();
};

#endif // COMPRESSIBLESWALLLSOBC_H
