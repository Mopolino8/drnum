#ifndef CARTESIANSTANDARDFLUXITERATOR_H
#define CARTESIANSTANDARDFLUXITERATOR_H

#include "patchiterator.h"
#include "cartesianpatchoperation.h"

class CartesianStandardIterator : public PatchIterator
{

private: // attributes

  CartesianPatchOperation* m_PatchOperation;
  size_t m_I1, m_I2;
  size_t m_J1, m_J2;
  size_t m_K1, m_K2;

public: // methods

  CartesianStandardIterator(CartesianPatchOperation *patch_operation);
  virtual void compute(real factor);
  void setDomain(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);
  void setI1(size_t i1) { m_I1 = i1; }
  void setI2(size_t i2) { m_I2 = i2; }
  void setJ1(size_t j1) { m_J1 = j1; }
  void setJ2(size_t j2) { m_J2 = j2; }
  void setK1(size_t k1) { m_K1 = k1; }
  void setK2(size_t k2) { m_K2 = k2; }

};


inline CartesianStandardIterator::CartesianStandardIterator(CartesianPatchOperation *patch_operation)
  : PatchIterator(patch_operation->patch())
{
  m_PatchOperation = patch_operation;
  m_I1 = 0;
  m_J1 = 0;
  m_K1 = 0;
  m_I2 = patch_operation->patch()->sizeI();
  m_J2 = patch_operation->patch()->sizeJ();
  m_K2 = patch_operation->patch()->sizeK();
}

inline void CartesianStandardIterator::compute(real factor)
{
  m_PatchOperation->compute(factor, m_I1, m_J1, m_K1, m_I2, m_J2, m_K2);
}

inline void CartesianStandardIterator::setDomain(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  m_I1 = i1;
  m_J1 = j1;
  m_K1 = k1;
  m_I2 = i2;
  m_J2 = j2;
  m_K2 = k2;
}

#endif // CARTESIANSTANDARDFLUXITERATOR_H
