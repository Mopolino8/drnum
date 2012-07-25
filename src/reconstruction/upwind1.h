#ifndef UPWIND1_H
#define UPWIND1_H

#include "cartesianpatch.h"

struct Upwind1
{
  static real xp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real xm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real yp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real ym(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real zp(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
  static real zm(CartesianPatch *P, size_t i_f, size_t i_v, size_t i, size_t j, size_t k);
};

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
