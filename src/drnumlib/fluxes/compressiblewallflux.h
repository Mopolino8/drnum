#ifndef COMPRESSIBLEWALLFLUX_H
#define COMPRESSIBLEWALLFLUX_H

#include "cartesianpatch.h"

template <typename TReconstruction, typename TGas>
class CompressibleWallFlux
{

  TReconstruction m_Reconstruction;

public:

  template <typename PATCH> CUDA_DH void xWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i-1, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*p;
    flux[2] += A*mu*v*patch->idx();
    flux[3] += A*mu*w*patch->idx();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void yWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j-1, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*mu*u*patch->idy();
    flux[2] += A*p;
    flux[3] += A*mu*w*patch->idy();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void zWallP(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k-1, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*mu*u*patch->idz();
    flux[2] += A*mu*v*patch->idz();
    flux[3] += A*p;
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void xWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] += A*p;
    flux[2] -= A*mu*v*patch->idx();
    flux[3] -= A*mu*w*patch->idx();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void yWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] -= A*mu*u*patch->idy();
    flux[2] += A*p;
    flux[3] -= A*mu*w*patch->idy();
    countFlops(2);
  }

  template <typename PATCH> CUDA_DH void zWallM(PATCH *patch, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    real var[5];
    patch->getVar(0, i, j, k, var);
    COMPR_VARS;
    real mu = TGas::mu(var);
    flux[1] -= A*mu*u*patch->idz();
    flux[2] -= A*mu*v*patch->idz();
    flux[3] += A*p;
    countFlops(2);
  }

};

#endif // COMPRESSIBLEWALLFLUX_H
