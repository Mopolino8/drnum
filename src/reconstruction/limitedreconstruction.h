#ifndef LIMITED_RECONSTRUCTION_H
#define LIMITED_RECONSTRUCTION_H

#include "cartesianpatch.h"

struct LimitedReconstruction
{
  static const real m_Epsilon;
  static real rx(CartesianPatch *P, size_t i_f, size_t i_v, size_t i1, size_t i2, size_t j, size_t k);
  static real ry(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j1, size_t j2, size_t k);
  static real rz(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k1, size_t k2);
};

inline real LimitedReconstruction::rx(CartesianPatch *P, size_t i_f, size_t i_v, size_t i1, size_t i2, size_t j, size_t k)
{
  countFlops(3);
  return CHECKED_REAL(   (P->f(i_f, i_v, i1, j, k) - P->f(i_f, i_v, 2*i1 - i2, j, k))
                       / nonZero(P->f(i_f, i_v, i2, j, k) - P->f(i_f, i_v, i1, j, k), m_Epsilon));
}

inline real LimitedReconstruction::ry(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j1, size_t j2, size_t k)
{
  countFlops(3);
  return CHECKED_REAL((P->f(i_f, i_v, i, j1, k) - P->f(i_f, i_v, i, 2*j1 - j2, k))/nonZero(P->f(i_f, i_v, i, j2, k) - P->f(i_f, i_v, i, j1, k), m_Epsilon));
}

inline real LimitedReconstruction::rz(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k1, size_t k2)
{
  countFlops(3);
  return CHECKED_REAL((P->f(i_f, i_v, i, j, k1) - P->f(i_f, i_v, i, j, 2*k1 - k2))/nonZero(P->f(i_f, i_v, i, j, k2) - P->f(i_f, i_v, i, j, k1), m_Epsilon));
}


struct FirstOrder  { static real lim(real)   { return 0; } };
struct SecondOrder { static real lim(real)   { return 1; } };
struct MinMod      { static real lim(real r) { countFlops(0); return CHECKED_REAL(max(0.0, min(1.0, r))); } };
struct SuperBee    { static real lim(real r) { countFlops(1); return CHECKED_REAL(max(0.0, max(min(2.0*r, 1.0), min(r, 2.0)))); } };
struct VanAlbada1  { static real lim(real r) { countFlops(5); return CHECKED_REAL((r*r + r)/(r*r + 1)); } };
struct VanAlbada2  { static real lim(real r) { countFlops(4); return CHECKED_REAL(2*r/(r*r + 1)); } };
struct VanLeer     { static real lim(real r) { countFlops(3); return CHECKED_REAL((r + fabs(r))/(1 + fabs(r))); } };

#endif // LIMITED_RECONSTRUCTION_H
