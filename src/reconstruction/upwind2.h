#ifndef UPWIND2_H
#define UPWIND2_H

#include "cartesianpatch.h"
#include "limitedreconstruction.h"

template <class TLimiter>
struct Upwind2 : public LimitedReconstruction
{
  static real xp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real xm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real yp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real ym(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real zp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real zm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
};


template <class TLimiter>
inline real Upwind2<TLimiter>::xp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (i > 1) {
    return P->f(var, i-1, j, k) + 0.5*TLimiter::lim(rx(P, var, i-1, i, j, k))*(P->f(var, i-1, j, k) - P->f(var, i-2, j, k));
  } else {
    return P->f(var, i-1, j, k) + 0.5*TLimiter::lim(rx(P, var, i, i+1, j, k))*(P->f(var, i, j, k) - P->f(var, i-1, j, k));
  }
}

template <class TLimiter>
inline real Upwind2<TLimiter>::xm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (i < P->sizeI() - 2) {
    return P->f(var, i, j, k) + 0.5*TLimiter::lim(rx(P, var, i, i-1, j, k))*(P->f(var, i, j, k) - P->f(var, i+1, j, k));
  } else {
    return P->f(var, i, j, k) + 0.5*TLimiter::lim(rx(P, var, i-1, i-2, j, k))*(P->f(var, i, j, k) - P->f(var, i-1, j, k));
  }
}


template <class TLimiter>
inline real Upwind2<TLimiter>::yp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (j > 1) {
    return P->f(var, i, j-1, k) + 0.5*TLimiter::lim(rx(P, var, i, j-1, j, k))*(P->f(var, i, j-1, k) - P->f(var, i, j-2, k));
  } else {
    return P->f(var, i, j-1, k) + 0.5*TLimiter::lim(rx(P, var, i, j, j+1, k))*(P->f(var, i, j, k) - P->f(var, i, j-1, k));
  }
}

template <class TLimiter>
inline real Upwind2<TLimiter>::ym(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (j < P->sizeJ() - 2) {
    return P->f(var, i, j, k) + 0.5*TLimiter::lim(rx(P, var, i, j, j-1, k))*(P->f(var, i, j, k) - P->f(var, i, j+1, k));
  } else {
    return P->f(var, i, j, k) + 0.5*TLimiter::lim(rx(P, var, i, j-1, j-2, k))*(P->f(var, i, j, k) - P->f(var, i, j-1, k));
  }
}


template <class TLimiter>
inline real Upwind2<TLimiter>::zp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (k > 1) {
    return P->f(var, i, j, k-1) + 0.5*TLimiter::lim(rx(P, var, i, j, k-1, k))*(P->f(var, i, j, k-1) - P->f(var, i, j, k-2));
  } else {
    return P->f(var, i, j, k-1) + 0.5*TLimiter::lim(rx(P, var, i, j, k, k+1))*(P->f(var, i, j, k) - P->f(var, i, j, k-1));
  }
}

template <class TLimiter>
inline real Upwind2<TLimiter>::zm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (k < P->sizeK() - 2) {
    return P->f(var, i, j, k) + 0.5*TLimiter::lim(rx(P, var, i, j, k, k-1))*(P->f(var, i, j, k) - P->f(var, i, j, k+1));
  } else {
    return P->f(var, i, j, k) + 0.5*TLimiter::lim(rx(P, var, i, j, k-1, k-2))*(P->f(var, i, j, k) - P->f(var, i, j, k-1));
  }
}

#endif // UPWIND2_H
