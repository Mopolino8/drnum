#ifndef COMPRESSIBLECARTESIANPATCH_H
#define COMPRESSIBLECARTESIANPATCH_H

#include "cartesianpatch.h"
#include "compressibleobject.h"

class CompressibleCartesianPatch : public CartesianPatch, public CompressibleObject
{

public: // methods

  CompressibleCartesianPatch();

  void setState(size_t i, size_t j, size_t k, real p, real T);
  void getState(size_t i, size_t j, size_t k, real &p, real &T);
  void getState(size_t i, size_t j, size_t k, real &p, real &u, real &v, real& w, real &T);

};

#endif // COMPRESSIBLECARTESIANPATCH_H
