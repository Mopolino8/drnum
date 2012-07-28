#ifndef UPWIND1_H
#define UPWIND1_H

#include "cartesianpatch.h"

struct Upwind1
{
  static void project(CartesianPatch *patch, real* var, size_t i_field, size_t num_vars,
                      size_t i1, size_t j1, size_t k1,
                      size_t,    size_t,    size_t);

  static real xp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real xm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real yp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real ym(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real zp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real zm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
};

inline void Upwind1::project(CartesianPatch *patch, real *var, size_t i_field, size_t num_vars,
                             size_t i1, size_t j1, size_t k1,
                             size_t   , size_t   , size_t)
{
  for (size_t i_var = 0; i_var < num_vars; ++i_var) {
    var[i_var] = patch->f(i_field, i_var, i1, j1, k1);
  }
}

inline real Upwind1::xp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  return P->f(i_f, i_v, i-1, j, k);
}

inline real Upwind1::xm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  return P->f(i_f, i_v, i, j, k);
}

inline real Upwind1::yp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  return P->f(i_f, i_v, i, j-1, k);
}

inline real Upwind1::ym(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  return P->f(i_f, i_v, i, j, k);
}

inline real Upwind1::zp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  return P->f(i_f, i_v, i, j, k-1);
}

inline real Upwind1::zm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k)
{
  return P->f(i_f, i_v, i, j, k);
}

#endif // UPWIND1_H
