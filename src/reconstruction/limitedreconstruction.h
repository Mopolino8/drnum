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

  static real r(CartesianPatch *patch, size_t i_f, size_t i_v, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);
};


inline real LimitedReconstruction::r(CartesianPatch *patch, size_t i_f, size_t i_v, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  countFlops(3);
  REGREAL delta1 = patch->f(i_f, i_v, i1, j1, k1) - patch->f(i_f, i_v, 2*i1-i2, 2*j1-j2, 2*k1-k2);
  REGREAL delta2 = patch->f(i_f, i_v, i2, j2, k2) - patch->f(i_f, i_v, i1, j1, k1);
  return CHECKED_REAL(delta1/nonZero(delta2, m_Epsilon));
}


struct FirstOrder  { static real lim(real)   { return 0; } };
struct SecondOrder { static real lim(real)   { return 1; } };
struct MinMod      { static real lim(real r) { countFlops(0); return CHECKED_REAL(max(real(0), min(real(1), r))); } };
struct SuperBee    { static real lim(real r) { countFlops(1); return CHECKED_REAL(max(real(0), max(min(real(2)*r, real(1)), min(r, real(2))))); } };
struct VanAlbada1  { static real lim(real r) { countFlops(5); return CHECKED_REAL((r*r + r)/(r*r + 1)); } };
struct VanAlbada2  { static real lim(real r) { countFlops(4); return CHECKED_REAL(2*r/(r*r + 1)); } };
struct VanLeer     { static real lim(real r) { countFlops(3); return CHECKED_REAL((r + fabs(r))/(1 + fabs(r))); } };

#endif // LIMITED_RECONSTRUCTION_H
