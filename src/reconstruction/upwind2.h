#ifndef UPWIND2_H
#define UPWIND2_H

#include "cartesianpatch.h"
#include "limitedreconstruction.h"

template <class TLimiter>
struct Upwind2 : public LimitedReconstruction
{
  void project(CartesianPatch *patch, real* var, size_t i_field, size_t num_vars,
               size_t i1, size_t j1, size_t k1,
               size_t i2, size_t j2, size_t k2,
               real = 0, real = 0, real = 0,
               real = 0, real = 0, real = 0);
};


template <class TLimiter>
inline void Upwind2<TLimiter>::project(CartesianPatch *patch, real *var, size_t i_field, size_t num_vars,
                                       size_t i1, size_t j1, size_t k1,
                                       size_t i2, size_t j2, size_t k2,
                                       real, real, real,
                                       real, real, real)
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
      //real lim = TLimiter::lim(r(patch, i_field, i_var, i1, j1, k1, i2, j2, k2));
      //var[i_var] += 0.5*lim*(patch->f(i_field, i_var, i1, j1, k1) - patch->f(i_field, i_var, i0, j0, k0));
      real delta01 = patch->f(i_field, i_var, i1, j1, k1) - patch->f(i_field, i_var, i0, j0, k0);
      real delta12 = patch->f(i_field, i_var, i2, j2, k2) - patch->f(i_field, i_var, i1, j1, k1);
      real delta = 0;
      if (delta01*delta12 > 0) {
        if (fabs(delta01) > fabs(delta12)) {
          delta = delta12;
        } else {
          delta = delta01;
        }
      }
      var[i_var] += 0.5*delta;
    }
  }
}

#endif // UPWIND2_H
