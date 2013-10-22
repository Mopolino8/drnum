#ifndef COMPRESSIBLESIMPLESWALLLSOBC_H
#define COMPRESSIBLESIMPLESWALLLSOBC_H

class CompressibleSimpleSWallLSOBC;

#include "levelsetobjectbc.h"

class CompressibleSimpleSWallLSOBC : public LevelSetObjectBC
{
public:

    CompressibleSimpleSWallLSOBC(size_t field,
                                 LevelSetObject* levelset_object,
                                 size_t abuse_field);

    virtual void operator()();
};

#endif // COMPRESSIBLESIMPLESWALLLSOBC_H
