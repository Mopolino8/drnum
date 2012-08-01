#ifndef CARTESIANSTANDARDFLUXITERATOR_H
#define CARTESIANSTANDARDFLUXITERATOR_H

#include "patchiterator.h"
#include "cartesianpatchoperation.h"

class CartesianStandardIterator : public PatchIterator
{

private: // attributes

  CartesianPatchOperation* m_PatchOperation;

public: // methods

  CartesianStandardIterator(CartesianPatchOperation *patch_operation) : PatchIterator(patch_operation->patch()) { m_PatchOperation = patch_operation; }
  virtual void compute(real factor);

};


inline void CartesianStandardIterator::compute(real factor)
{
  transformShapes();
  size_t i2 = m_PatchOperation->patch()->sizeI();
  size_t j2 = m_PatchOperation->patch()->sizeJ();
  size_t k2 = m_PatchOperation->patch()->sizeK();
  m_PatchOperation->compute(factor, 0, 0, 0, i2, j2, k2);
}

#endif // CARTESIANSTANDARDFLUXITERATOR_H
