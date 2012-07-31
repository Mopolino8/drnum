#ifndef CARTESIANDIRECTIONALPATCHOPERATION_H
#define CARTESIANDIRECTIONALPATCHOPERATION_H

#include "cartesianpatchoperation.h"

template <unsigned int DIM, class TFlux>
class CartesianDirectionalPatchOperation : public CartesianPatchOperation
{

  TFlux *m_Flux;

public:

  CartesianDirectionalPatchOperation(CartesianPatch *patch, TFlux* flux);

  virtual void compute(real factor, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);

};


template <unsigned int DIM, class TFlux>
CartesianDirectionalPatchOperation<DIM, TFlux>::CartesianDirectionalPatchOperation(CartesianPatch *patch, TFlux* flux) : CartesianPatchOperation(patch)
{
  m_Flux = flux;
}


template <unsigned int DIM, class TFlux>
void CartesianDirectionalPatchOperation<DIM, TFlux>::compute(real factor, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  checkResFieldSize(i1, j1, k1, i2, j2, k2);

  real Ax = patch()->dy()*patch()->dz();
  real Ay = patch()->dx()*patch()->dz();
  real Az = patch()->dx()*patch()->dy();

  real flux[5];

  for (size_t i_res = 0; i_res < DIM*m_ResLength; ++i_res) {
    m_Res[i_res] = 0;
  }
  countFlops(3);

  // compute main block
  //
  // .. x direction
  for (size_t j = 0; j < patch()->sizeJ(); ++j) {
    for (size_t k = 0; k < patch()->sizeK(); ++k) {
      for (size_t i = 0; i < patch()->sizeI(); ++i) {
        if (i > 0) {
          fill(flux, 5, 0);
          m_Flux->x(patch(), i, j, k, Ax, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i-1, j, k)] -= flux[i_var];
            m_Res[resIndex(i_var, i, j, k)]   += flux[i_var];
          }
          countFlops(2*DIM);
        }
      }
    }
  }

  // .. y direction
  for (size_t i = 0; i < patch()->sizeI(); ++i) {
    for (size_t k = 0; k < patch()->sizeK(); ++k) {
      for (size_t j = 0; j < patch()->sizeJ(); ++j) {
        if (j > 0) {
          fill(flux, 5, 0);
          m_Flux->y(patch(), i, j, k, Ay, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i, j-1, k)] -= flux[i_var];
            m_Res[resIndex(i_var, i, j, k)]   += flux[i_var];
          }
          countFlops(2*DIM);
        }
      }
    }
  }

  // .. z direction
  for (size_t i = 0; i < patch()->sizeI(); ++i) {
    for (size_t j = 0; j < patch()->sizeJ(); ++j) {
      for (size_t k = 0; k < patch()->sizeK(); ++k) {
        // z direction
        if (k > 0) {
          fill(flux, 5, 0);
          m_Flux->z(patch(), i, j, k, Az, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i, j, k-1)] -= flux[i_var];
            m_Res[resIndex(i_var, i, j, k)]   += flux[i_var];
          }
          countFlops(2*DIM);
        }
      }
    }
  }

  // compute x walls
  //
  // .. left wall
  if (i1 == 0) {
    for (size_t j = j1; j < j2; ++j) {
      for (size_t k = k1; k < k2; ++k) {
        fill(flux, 5, 0);
        m_Flux->xWallM(patch(), i1, j, k, Ax, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i1, j, k)] += flux[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // .. right wall
  if (i2 == patch()->sizeI()) {
    for (size_t j = j1; j < j2; ++j) {
      for (size_t k = k1; k < k2; ++k) {
        fill(flux, 5, 0);
        m_Flux->xWallP(patch(), i2, j, k, Ax, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i2-1, j, k)] -= flux[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // compute y walls
  //
  // .. front wall
  if (j1 == 0) {
    for (size_t i = i1; i < i2; ++i) {
      for (size_t k = k1; k < k2; ++k) {
        fill(flux, 5, 0);
        m_Flux->yWallM(patch(), i, j1, k, Ay, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j1, k)] += flux[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // .. back wall
  if (j2 == patch()->sizeJ()) {
    for (size_t i = i1; i < i2; ++i) {
      for (size_t k = k1; k < k2; ++k) {
        fill(flux, 5, 0);
        m_Flux->yWallP(patch(), i, j2, k, Ay, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j2-1, k)] -= flux[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // compute z walls
  //
  // .. bottom wall
  if (k1 == 0) {
    for (size_t i = i1; i < i2; ++i) {
      for (size_t j = j1; j < j2; ++j) {
        fill(flux, 5, 0);
        m_Flux->zWallM(patch(), i, j, k1, Az, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j, k1)] += flux[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // .. top wall
  if (k2 == patch()->sizeK()) {
    for (size_t i = i1; i < i2; ++i) {
      for (size_t j = j1; j < j2; ++j) {
        fill(flux, 5, 0);
        m_Flux->zWallP(patch(), i, j, k2, Az, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j, k2-1)] -= flux[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // advance to next iteration level (time)
  factor /= patch()->dV();
  for (size_t i = 0; i < patch()->sizeI(); ++i) {
    for (size_t j = 0; j < patch()->sizeJ(); ++j) {
      for (size_t k = 0; k < patch()->sizeK(); ++k) {
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          patch()->f(0, i_var, i, j, k) = patch()->f(1, i_var, i, j, k) + factor*m_Res[resIndex(i_var, i, j, k)];
        }
      }
    }
  }
}

#endif // CARTESIANDIRECTIONALPATCHOPERATION_H
