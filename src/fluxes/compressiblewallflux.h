#ifndef COMPRESSIBLEWALLFLUX_H
#define COMPRESSIBLEWALLFLUX_H

#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class CompressibleWallFlux
{

  TReconstruction* m_Reconstruction;

public:

  CompressibleWallFlux(TReconstruction* reconstruction) { m_Reconstruction = reconstruction; }

  void xWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void yWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void zWallP(CartesianPatch *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void xWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void yWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void zWallM(CartesianPatch *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

};

#define COMPR_VARS \
  real r  = var[0]; \
  real ir = CHECKED_REAL(1.0/r); \
  real ru = var[1]; \
  real rv = var[2]; \
  real rw = var[3]; \
  real u  = ru*ir; \
  real v  = rv*ir; \
  real w  = rw*ir; \
  real rE = var[4]; \
  real T  = CHECKED_REAL((rE*ir - 0.5*(u*u + v*v + w*w))/TGas::cv(var)); \
  real p  = r*TGas::R(var)*T; \
  countFlops(15);

template <typename TReconstruction, typename TGas>
inline void CompressibleWallFlux<TReconstruction, TGas>::xWallP(CartesianPatch *patch,
                                                                size_t i, size_t j, size_t k,
                                                                real x, real y, real z,
                                                                real A, real* flux)
{
  real var[5];
  m_Reconstruction->project(patch, var, 0, 5, i-1, j, k, i, j, k, x - patch->dx(), y, z, x, y, z);
  COMPR_VARS;
  flux[1] += A*p;
  countFlops(2);
}

template <typename TReconstruction, typename TGas>
inline void CompressibleWallFlux<TReconstruction, TGas>::xWallM(CartesianPatch *patch,
                                                                size_t i, size_t j, size_t k,
                                                                real x, real y, real z,
                                                                real A, real* flux)
{
  real var[5];
  m_Reconstruction->project(patch, var, 0, 5, i, j, k, i-1, j, k, x, y, z, x - patch->dx(), y, z);
  COMPR_VARS;
  flux[1] += A*p;
  countFlops(2);
}


template <typename TReconstruction, typename TGas>
inline void CompressibleWallFlux<TReconstruction, TGas>::yWallP(CartesianPatch *patch,
                                                                size_t i, size_t j, size_t k,
                                                                real x, real y, real z,
                                                                real A, real* flux)
{
  real var[5];
  m_Reconstruction->project(patch, var, 0, 5, i, j-1, k, i, j, k, x, y - patch->dy(), z, x, y, z);
  COMPR_VARS;
  flux[2] += A*p;
  countFlops(2);
}

template <typename TReconstruction, typename TGas>
inline void CompressibleWallFlux<TReconstruction, TGas>::yWallM(CartesianPatch *patch,
                                                                size_t i, size_t j, size_t k,
                                                                real x, real y, real z,
                                                                real A, real* flux)
{
  real var[5];
  m_Reconstruction->project(patch, var, 0, 5, i, j, k, i, j-1, k, x, y, z, x, y - patch->dy(), z);
  COMPR_VARS;
  flux[2] += A*p;
  countFlops(2);
}


template <typename TReconstruction, typename TGas>
inline void CompressibleWallFlux<TReconstruction, TGas>::zWallP(CartesianPatch *patch,
                                                                size_t i, size_t j, size_t k,
                                                                real x, real y, real z,
                                                                real A, real* flux)
{
  real var[5];
  m_Reconstruction->project(patch, var, 0, 5, i, j, k-1, i, j, k, x, y, z - patch->dz(), x, y, z);
  COMPR_VARS;
  flux[3] += A*p;
  countFlops(2);
}

template <typename TReconstruction, typename TGas>
inline void CompressibleWallFlux<TReconstruction, TGas>::zWallM(CartesianPatch *patch,
                                                                size_t i, size_t j, size_t k,
                                                                real x, real y, real z,
                                                                real A, real* flux)
{
  real var[5];
  m_Reconstruction->project(patch, var, 0, 5, i, j, k, i, j, k-1, x, y, z, x, y, z - patch->dz());
  COMPR_VARS;
  flux[3] += A*p;
  countFlops(2);
}

#endif // COMPRESSIBLEWALLFLUX_H
