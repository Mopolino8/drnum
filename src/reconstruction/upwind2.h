#ifndef UPWIND2_H
#define UPWIND2_H

#include "cartesianpatch.h"
#include "limited_reconstruction.h"

template <class TLimiter, real TEpsilon>
struct Upwind2 : public LimitedReconstruction<TEpsilon>
{
  static real xp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real xm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real yp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real ym(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real zp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
  static real zm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k);
};

inline real Upwind2::xp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (i > 1) {
    return P->f(var, i-1, j, k) + 0.5*TLimiter::lim(rx(var, i-1, i, j, k))*(P->f(var, i-1, j, k) - P->f(var, i-2, j, k));
  } else {
    return P->f(var, i-1, j, k) + 0.5*TLimiter::lim(rx(var, i, i+1, j, k))*(P->f(var, i, j, k) - P->f(var, i-1, j, k));
  }
}

inline real Upwind2::xm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  if (i < P->sizeI() - 2) {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, i-1, j, k))*(f(var, i, j, k) - f(var, i+1, j, k));
  } else {


     //hier stimmts noch nicht !!!
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i+1, i, j, k))*(f(var, i+1, j, k) - f(var, i, j, k));
  }
}

inline real Upwind1::yp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  return P->f(var, i, j-1, k);
}

inline real Upwind1::ym(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  return P->f(var, i, j, k);
}

inline real Upwind1::zp(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  return P->f(var, i, j, k-1);
}

inline real Upwind1::zm(CartesianPatch *P, real *var, size_t i, size_t j, size_t k)
{
  return P->f(var, i, j, k);
}


inline real CartesianPatch::projUxP(real *var, size_t i, size_t j, size_t k)
{
  if (i > 0) {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i-1, i, j, k))*(f(var, i, j, k) - f(var, i-1, j, k));
  } else {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, i+1, j, k))*(f(var, i+1, j, k) - f(var, i, j, k));
  }
}

inline real CartesianPatch::projUxM(real *var, size_t i, size_t j, size_t k)
{
  if (i < m_NumI - 1) {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i+1, i, j, k))*(f(var, i, j, k) - f(var, i+1, j, k));
  } else {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, i-1, j, k))*(f(var, i-1, j, k) - f(var, i, j, k));
  }
}


inline real CartesianPatch::projUyP(real *var, size_t i, size_t j, size_t k)
{
  if (j > 0) {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j-1, j, k))*(f(var, i, j, k) - f(var, i, j-1, k));
  } else {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j, j+1, k))*(f(var, i, j+1, k) - f(var, i, j, k));
  }
}

inline real CartesianPatch::projUyM(real *var, size_t i, size_t j, size_t k)
{
  if (j < m_NumJ - 1) {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j+1, j, k))*(f(var, i, j, k) - f(var, i, j+1, k));
  } else {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j, j-1, k))*(f(var, i, j-1, k) - f(var, i, j, k));
  }
}


inline real CartesianPatch::projUzP(real *var, size_t i, size_t j, size_t k)
{
  if (k > 0) {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j, k-1, k))*(f(var, i, j, k) - f(var, i, j, k-1));
  } else {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j, k, k+1))*(f(var, i, j, k+1) - f(var, i, j, k));
  }
}

inline real CartesianPatch::projUzM(real *var, size_t i, size_t j, size_t k)
{
  if (i < m_NumI - 1) {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j, k+1, j))*(f(var, i, j, k) - f(var, i, j, k+1));
  } else {
    return f(var, i, j, k) + 0.5*TLimiter::lim(rx(var, i, j, k, k-1))*(f(var, i, j, k-1) - f(var, i, j, k));
  }
}
*/

#endif // UPWIND2_H
