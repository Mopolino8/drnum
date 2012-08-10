#ifndef CARTESIANPATCHOPERATION_H
#define CARTESIANPATCHOPERATION_H

#include "cartesianpatch.h"

class CartesianPatchOperation
{

private: // attributes

  CartesianPatch *m_Patch;

protected: // attributes

  real*  m_Res;
  size_t m_ResLength;
  size_t m_I1;
  size_t m_J1;
  size_t m_K1;
  size_t m_SizeI;
  size_t m_SizeJ;
  size_t m_SizeK;


protected: // methods

  void checkResFieldSize(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);


public: // methods

  CartesianPatchOperation(CartesianPatch *P);
  virtual ~CartesianPatchOperation();

  CartesianPatch* patch() { return m_Patch; }

  size_t resIndex(size_t i_var, size_t i, size_t j, size_t k) { return m_ResLength*i_var + (i-m_I1)*m_SizeJ*m_SizeK + (j-m_J1)*m_SizeK + (k-m_K1); }

  virtual void compute(real factor, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2) = 0;

};


inline CartesianPatchOperation::CartesianPatchOperation(CartesianPatch *P)
{
  m_Patch = P;
  m_Res = NULL;
  m_ResLength = 0;
}

inline CartesianPatchOperation::~CartesianPatchOperation()
{
  delete [] m_Res;
}

inline void CartesianPatchOperation::checkResFieldSize(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  m_SizeI = i2 - i1;
  m_SizeJ = j2 - j1;
  m_SizeK = k2 - k1;
  m_I1 = i1;
  m_J1 = j1;
  m_K1 = k1;
  size_t new_length = m_SizeI*m_SizeJ*m_SizeK;
  if (new_length > m_ResLength) {
    delete [] m_Res;
    m_Res = new real [new_length*patch()->numVariables()];
    m_ResLength = new_length;
  }
}



#endif // CARTESIANPATCHOPERATION_H
