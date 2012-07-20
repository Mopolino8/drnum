#ifndef LIMITED_RECONSTRUCTION_H
#define LIMITED_RECONSTRUCTION_H

#include "cartesianpatch.h"

template <real TEpsilon>
class LimitedReconstruction
{

  real rx(CartesianPatch *P, real *var, size_t i1, size_t i2, size_t j, size_t k);
  real ry(CartesianPatch *P, real *var, size_t i, size_t j1, size_t j2, size_t k);
  real rz(CartesianPatch *P, real *var, size_t i, size_t j, size_t k1, size_t k2);

};

inline LimitedReconstruction::rx(CartesianPatch *P, real *var, size_t i1, size_t i2, size_t j, size_t k)
{
  return (P->f(var, i1, j, k) - P->f(var, 2*i1 - i2, j, k))/nonZero(P->f(var, i2, j, k) - P->f(var, i1, j, k), TEpsilon);
}

inline real LimitedReconstruction::ry(CartesianPatch *P, real *var, size_t i, size_t j1, size_t j2, size_t k)
{
  return (P->f(var, i, j1, k) - P->f(var, i, 2*j1 - j2, k))/nonZero(P->f(var, i, j2, k) - P->f(var, i, j1, k), TEpsilon);
}

inline real LimitedReconstruction::rz(CartesianPatch *P, real *var, size_t i, size_t j, size_t k1, size_t k2)
{
  return (P->f(var, i, j, k1) - P->f(var, i, j, 2*k1 - k2))/nonZero(P->f(var, i, j, k2) - P->f(var, i, j, k1), TEpsilon);
}

#endif // LIMITED_RECONSTRUCTION_H
