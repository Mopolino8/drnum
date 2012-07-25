#ifndef COMPRESSIBLECARTESIANPATCHOPERATION_H
#define COMPRESSIBLECARTESIANPATCHOPERATION_H

#include "cartesianpatchoperation.h"
#include "fluxes/compressibleflux.h"

template <class TReconstruction>
class CompressibleCartesianPatchOperation : public CartesianPatchOperation
{

private: // attributes

  CompressibleFlux m_Flux;

protected: // methods

  void setPointers(CompressibleFlux *flux);
  void computeWalls(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);

};

template <class TReconstruction>
inline void CompressibleCartesianPatchOperation<TReconstruction>::setPointers(CompressibleFlux *flux)
{
  flux->setRho (getVariable(i_new, 0));
  flux->setRhou(getVariable(i_new, 1));
  flux->setRhov(getVariable(i_new, 2));
  flux->setRhow(getVariable(i_new, 3));
  flux->setRhoE(getVariable(i_new, 4));
  flux->setResRho (getVariable(i_res, 0));
  flux->setResRhou(getVariable(i_res, 1));
  flux->setResRhov(getVariable(i_res, 2));
  flux->setResRhow(getVariable(i_res, 3));
  flux->setResRhoE(getVariable(i_res, 4));
}

template <class TReconstruction>
inline void CompressibleCartesianPatchOperation<TReconstruction>::computeWalls(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  if (i1 == 0) {
    for (size_t j = j1; j < j2; ++j) {
      for (size_t k = k1; k < k2; ++k) {
        m_Flux.xWallM(patch(), i, j, k, A, flux);
      }
    }
  }
}

#endif // COMPRESSIBLECARTESIANPATCHOPERATION_H
