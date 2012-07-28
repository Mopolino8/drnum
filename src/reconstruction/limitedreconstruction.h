#ifndef LIMITED_RECONSTRUCTION_H
#define LIMITED_RECONSTRUCTION_H

#include "cartesianpatch.h"

/**
 * Base class for reconstruction methods.
 * The implementation is probably not optimal at the moment,
 * especially the unified methods (for x, y, z directions) are most likely
 * not very efficient due to too many if statements. Later on this can be split again,
 * but for now it should reduce the chance of errors significantly.
 * @todo look at efficiency
 */
struct LimitedReconstruction
{
  static const real m_Epsilon;

  /**
   * Check if an (i,j,k) is inside a Cartesian patch.
   * Attention only the upper limit will be checked (unsigned data type).
   * @param patch The patch to check
   * @param i first index
   * @param j second index
   * @param k third index
   * @return true if it is a valid (i,j,k) triple, false otherwise
   */
  static bool checkRange(CartesianPatch *patch, size_t i, size_t j, size_t k);

  static real r(CartesianPatch *patch, size_t i_f, size_t i_v, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);

  static real rx(CartesianPatch *P, size_t i_f, size_t i_v, size_t i1, size_t i2, size_t j, size_t k);
  static real ry(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j1, size_t j2, size_t k);
  static real rz(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k1, size_t k2);
};


inline bool LimitedReconstruction::checkRange(CartesianPatch *patch, size_t i, size_t j, size_t k)
{
  if (i >= patch->sizeI() || j >= patch->sizeJ() || k >= patch->sizeK()) {
    return false;
  }
  return true;
}

inline real LimitedReconstruction::r(CartesianPatch *patch, size_t i_f, size_t i_v, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  countFlops(3);
  REGREAL delta1 = patch->f(i_f, i_v, i1, j1, k1) - patch->f(i_f, i_v, 2*i1-i2, 2*j1-j2, 2*k1-k2);
  REGREAL delta2 = patch->f(i_f, i_v, i2, j2, k2) - patch->f(i_f, i_v, i1, j1, k1);
  return CHECKED_REAL(delta1/nonZero(delta2, m_Epsilon));
}


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
struct MinMod      { static real lim(real r) { countFlops(0); return CHECKED_REAL(max(real(0), min(real(1), r))); } };
struct SuperBee    { static real lim(real r) { countFlops(1); return CHECKED_REAL(max(real(0), max(min(real(2)*r, real(1)), min(r, real(2))))); } };
struct VanAlbada1  { static real lim(real r) { countFlops(5); return CHECKED_REAL((r*r + r)/(r*r + 1)); } };
struct VanAlbada2  { static real lim(real r) { countFlops(4); return CHECKED_REAL(2*r/(r*r + 1)); } };
struct VanLeer     { static real lim(real r) { countFlops(3); return CHECKED_REAL((r + fabs(r))/(1 + fabs(r))); } };

#endif // LIMITED_RECONSTRUCTION_H
