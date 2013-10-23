#ifndef COMPRESSIBLEEULERLSOBC_H
#define COMPRESSIBLEEULERLSOBC_H

class CompressibleEulerLSOBC;

#include "levelsetobjectbc.h"

class CompressibleEulerLSOBC : public LevelSetObjectBC
{
public:

  CompressibleEulerLSOBC(size_t field,
                         LevelSetObject* levelset_object,
                         size_t abuse_field);

  virtual void operator()();
};

#endif // COMPRESSIBLEEULERLSOBC_H
