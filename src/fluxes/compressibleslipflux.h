#ifndef COMPRESSIBLEWALLFLUX_H
#define COMPRESSIBLEWALLFLUX_H

#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class CompressibleSlipFlux
{

  TReconstruction m_Reconstruction;

public:

  template <typename PATCH> CUDA_DH void xWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    //m_Reconstruction->project(patch, var, 0, 5, i-1, j, k, i, j, k, x - patch->dx(), y, z, x, y, z);
    patch->getVar(0, i-1, j, k, var);
    COMPR_VARS;
    flux[1] += A*p;
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void yWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    //m_Reconstruction->project(patch, var, 0, 5, i, j-1, k, i, j, k, x, y - patch->dy(), z, x, y, z);
    patch->getVar(0, i, j-1, k, var);
    COMPR_VARS;
    flux[2] += A*p;
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void zWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    //m_Reconstruction->project(patch, var, 0, 5, i, j, k-1, i, j, k, x, y, z - patch->dz(), x, y, z);
    patch->getVar(0, i, j, k-1, var);
    COMPR_VARS;
    flux[3] += A*p;
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void xWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    //m_Reconstruction->project(patch, var, 0, 5, i, j, k, i-1, j, k, x, y, z, x - patch->dx(), y, z);
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    flux[1] += A*p;
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void yWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    //m_Reconstruction->project(patch, var, 0, 5, i, j, k, i, j-1, k, x, y, z, x, y - patch->dy(), z);
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    flux[2] += A*p;
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void zWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    //m_Reconstruction->project(patch, var, 0, 5, i, j, k, i, j, k-1, x, y, z, x, y, z - patch->dz());
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    flux[3] += A*p;
    countFlops(2);
  }

};

#endif // COMPRESSIBLEWALLFLUX_H
