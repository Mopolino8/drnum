#ifndef UPWIND2_H
#define UPWIND2_H

#include "cartesianpatch.h"
#include "limitedreconstruction.h"

template <class TLimiter>
struct Upwind2 : public LimitedReconstruction
{
  static void project(CartesianPatch *patch, real* var, size_t i_field, size_t num_vars,
                      size_t i1, size_t j1, size_t k1,
                      size_t i2, size_t j2, size_t k2);

  static real xp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real xm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real yp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real ym(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real zp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real zm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
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
  if (checkRange(patch, i0, j0, k0) && checkRange(patch, i2, j2, k2)) {
    for (size_t i_var = 0; i_var < num_vars; ++i_var) {
      REGREAL lim = TLimiter::lim(r(patch, i_field, i_var, i1, j1, k1, i2, j2, k2));
      var[i_var] += 0.5*lim*(patch->f(i_field, i_var, i1, j1, k1) - patch->f(i_field, i_var, i0, j0, k0));
    }
  }
}


template <class TLimiter>
inline real Upwind2<TLimiter>::xp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  if (i > 1) {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i-1, j, k) + 0.5*TLimiter::lim(rx(P, i_f, i_v, i-1, i, j, k))*(P->f(i_f, i_v, i-1, j, k) - P->f(i_f, i_v, i-2, j, k)));
  } else {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i-1, j, k) + 0.5*TLimiter::lim(rx(P, i_f, i_v, i, i+1, j, k))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i-1, j, k)));
    //return P->f(i_f, i_v, i-1, j, k);
  }
}

template <class TLimiter>
inline real Upwind2<TLimiter>::xm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  if (i < P->sizeI() - 2) {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k) + 0.5*TLimiter::lim(rx(P, i_f, i_v, i, i-1, j, k))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i+1, j, k)));
    //if (i > 0) return CHECKED_REAL(P->f(i_f, i_v, i, j, k) + 0.5*TLimiter::lim(rx(P, i_f, i_v, i, i-1, j, k))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i+1, j, k)));
    //return P->f(i_f, i_v, i, j, k);
  } else {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k) + 0.5*TLimiter::lim(rx(P, i_f, i_v, i-1, i-2, j, k))*(P->f(i_f, i_v, i-1, j, k) - P->f(i_f, i_v, i, j, k)));
    //return P->f(i_f, i_v, i, j, k);
  }
}


template <class TLimiter>
inline real Upwind2<TLimiter>::yp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  if (j > 1) {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j-1, k) + 0.5*TLimiter::lim(ry(P, i_f, i_v, i, j-1, j, k))*(P->f(i_f, i_v, i, j-1, k) - P->f(i_f, i_v, i, j-2, k)));
  } else {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j-1, k) + 0.5*TLimiter::lim(ry(P, i_f, i_v, i, j, j+1, k))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i, j-1, k)));
    //return P->f(i_f, i_v, i, j-1, k);
  }
}

template <class TLimiter>
inline real Upwind2<TLimiter>::ym(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  if (j < P->sizeJ() - 2) {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k) + 0.5*TLimiter::lim(ry(P, i_f, i_v, i, j, j-1, k))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i, j+1, k)));
  } else {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k) + 0.5*TLimiter::lim(ry(P, i_f, i_v, i, j-1, j-2, k))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i, j-1, k)));
    //return P->f(i_f, i_v, i, j, k);
  }
}


template <class TLimiter>
inline real Upwind2<TLimiter>::zp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  if (k > 1) {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k-1) + 0.5*TLimiter::lim(rz(P, i_f, i_v, i, j, k-1, k))*(P->f(i_f, i_v, i, j, k-1) - P->f(i_f, i_v, i, j, k-2)));
  } else {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k-1) + 0.5*TLimiter::lim(rz(P, i_f, i_v, i, j, k, k+1))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i, j, k-1)));
    //return P->f(i_f, i_v, i, j, k-1);
  }
}

template <class TLimiter>
inline real Upwind2<TLimiter>::zm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  if (k < P->sizeK() - 2) {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k) + 0.5*TLimiter::lim(rz(P, i_f, i_v, i, j, k, k-1))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i, j, k+1)));
  } else {
    countFlops(4);
    return CHECKED_REAL(P->f(i_f, i_v, i, j, k) + 0.5*TLimiter::lim(rz(P, i_f, i_v, i, j, k-1, k-2))*(P->f(i_f, i_v, i, j, k) - P->f(i_f, i_v, i, j, k-1)));
    //return P->f(i_f, i_v, i, j, k);
  }
}

#endif // UPWIND2_H
