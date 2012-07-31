#ifndef UPWIND2_H
#define UPWIND2_H

#include "cartesianpatch.h"
#include "limitedreconstruction.h"

template <class TLimiter>
struct Upwind2 : public LimitedReconstruction
{
  void project(CartesianPatch *patch, real* var, size_t i_field, size_t num_vars,
               size_t i1, size_t j1, size_t k1,
               size_t i2, size_t j2, size_t k2);
};


template <class TLimiter>
inline void Upwind2<TLimiter>::project(CartesianPatch *patch, real *var, size_t i_field, size_t num_vars,
                                       size_t i1, size_t j1, size_t k1,
                                       size_t i2, size_t j2, size_t k2)
{
  /// @todo re-enable second order boundary reconstruction
  for (size_t i_var = 0; i_var < num_vars; ++i_var) {
    var[i_var] = patch->f(i_field, i_var, i1, j1, k1);
  }
  size_t i0 = 2*i1 - i2;
  size_t j0 = 2*j1 - j2;
  size_t k0 = 2*k1 - k2;
  if (patch->checkRange(i0, j0, k0) && patch->checkRange(i2, j2, k2)) {
    for (size_t i_var = 0; i_var < num_vars; ++i_var) {
      REGREAL lim = TLimiter::lim(r(patch, i_field, i_var, i1, j1, k1, i2, j2, k2));
      var[i_var] += 0.5*lim*(patch->f(i_field, i_var, i1, j1, k1) - patch->f(i_field, i_var, i0, j0, k0));
    }
  }
}

#endif // UPWIND2_H
