#ifndef COMPRESSIBLEEULERLSOBC_H
#define COMPRESSIBLEEULERLSOBC_H

class CompressibleEulerLSOBC;

#include "levelsetobjectbc.h"

class CompressibleEulerLSOBC : public LevelSetObjectBC
{
public:

  CompressibleEulerLSOBC(size_t field, LevelSetObject* levelset_object);

  virtual void operator()();
};

#endif // COMPRESSIBLEEULERLSOBC_H
