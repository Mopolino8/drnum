#ifndef UPWINDCENTRAL_H
#define UPWINDCENTRAL_H

#include "cartesianpatch.h"

template <class TLimiter>
struct UpwindCentral
{
  void project(CartesianPatch *patch, real* var, size_t i_field, size_t num_vars,
               size_t i1, size_t j1, size_t k1,
               size_t i2, size_t j2, size_t k2,
               real = 0, real = 0, real = 0,
               real = 0, real = 0, real = 0);
};


template <class TLimiter>
inline void UpwindCentral<TLimiter>::project(CartesianPatch *patch, real *var, size_t i_field, size_t num_vars,
                                             size_t i1, size_t j1, size_t k1,
                                             size_t i2, size_t j2, size_t k2,
                                             real, real, real,
                                             real, real, real)
{
  for (size_t i_var = 0; i_var < num_vars; ++i_var) {
    var[i_var] = patch->f(i_field, i_var, i1, j1, k1);
  }
  size_t i0 = 2*i1 - i2;
  size_t j0 = 2*j1 - j2;
  size_t k0 = 2*k1 - k2;
  if (patch->checkRange(i0, j0, k0) && patch->checkRange(i2, j2, k2)) {
    for (size_t i_var = 0; i_var < num_vars; ++i_var) {
      real delta01 = patch->f(i_field, i_var, i1, j1, k1) - patch->f(i_field, i_var, i0, j0, k0);
      real delta12 = patch->f(i_field, i_var, i2, j2, k2) - patch->f(i_field, i_var, i1, j1, k1);
      var[i_var] += 0.5*TLimiter::lim(delta01, delta12)*delta12;
    }
  }
}

#endif // UPWINDCENTRAL_H
