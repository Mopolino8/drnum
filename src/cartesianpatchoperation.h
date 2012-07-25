#ifndef CARTESIANPATCHOPERATION_H
#define CARTESIANPATCHOPERATION_H

#include "cartesianpatch.h"

class CartesianPatchOperation
{

private: // attributes

  CartesianPatch *m_Patch;

public: // methods

  CartesianPatchOperation(CartesianPatch *P) { m_Patch = P; }
  CartesianPatch* patch() { return m_Patch; }

  virtual void compute(real factor, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2) = 0;

};

#endif // CARTESIANPATCHOPERATION_H
