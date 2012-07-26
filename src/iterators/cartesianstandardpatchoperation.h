#ifndef CARTESIANSTANDARDPATCHOPERATION_H
#define CARTESIANSTANDARDPATCHOPERATION_H

#include "cartesianpatchoperation.h"

template <unsigned int DIM, class TFlux>
class CartesianStandardPatchOperation : public CartesianPatchOperation
{

private: // attributes

  TFlux  m_Flux;
  real*  m_Res;
  size_t m_ResLength;
  size_t m_I1;
  size_t m_J1;
  size_t m_K1;
  size_t m_SizeI;
  size_t m_SizeJ;
  size_t m_SizeK;


protected: // methods

  void checkResFieldSize(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);


public:

  CartesianStandardPatchOperation(CartesianPatch *patch);
  virtual ~CartesianStandardPatchOperation();

  size_t resIndex(size_t i_var, size_t i, size_t j, size_t k) { return m_ResLength*i_var + (i-m_I1)*m_SizeJ*m_SizeK + (j-m_J1)*m_SizeK + (k-m_K1); }
  virtual void compute(real factor, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2);

};


template <unsigned int DIM, class TFlux>
CartesianStandardPatchOperation<DIM, TFlux>::CartesianStandardPatchOperation(CartesianPatch *patch) : CartesianPatchOperation(patch)
{
  m_Res = NULL;
  m_ResLength = 0;
}


template <unsigned int DIM, class TFlux>
CartesianStandardPatchOperation<DIM, TFlux>::~CartesianStandardPatchOperation()
{
  delete [] m_Res;
}


template <unsigned int DIM, class TFlux>
void CartesianStandardPatchOperation<DIM, TFlux>::checkResFieldSize(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  m_SizeI = i2 - i1;
  m_SizeJ = j2 - j1;
  m_SizeK = k2 - k1;
  m_I1 = i1;
  m_J1 = j1;
  m_K1 = k1;
  size_t new_length = m_SizeI*m_SizeJ*m_SizeK;
  if (new_length > m_ResLength) {
    delete [] m_Res;
    m_Res = new real [new_length*patch()->numVariables()];
    m_ResLength = new_length;
  }
}


template <unsigned int DIM, class TFlux>
void CartesianStandardPatchOperation<DIM, TFlux>::compute(real factor, size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2)
{
  checkResFieldSize(i1, j1, k1, i2, j2, k2);

  real Ax = patch()->dy()*patch()->dz();
  real Ay = patch()->dx()*patch()->dz();
  real Az = patch()->dx()*patch()->dy();
  RealVec<DIM> flux;
  for (size_t i_res = 0; i_res < DIM*m_ResLength; ++i_res) {
    m_Res[i_res] = 0;
  }
  countFlops(3);

  // compute main block
  for (size_t i = 0; i < patch()->sizeI(); ++i) {
    for (size_t j = 0; j < patch()->sizeJ(); ++j) {
      for (size_t k = 0; k < patch()->sizeK(); ++k) {

        // x direction
        if (i > 0) {
          flux.fill(0);
          m_Flux.x(patch(), i, j, k, Ax, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            if (i > 0) {
              m_Res[resIndex(i_var, i-1, j, k)] -= flux.var[i_var];
              countFlops(1);
            }
            m_Res[resIndex(i_var, i, j, k)]   += flux.var[i_var];
          }
          countFlops(DIM);
        }

        // y direction
        if (j > 0) {
          flux.fill(0);
          m_Flux.y(patch(), i, j, k, Ay, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            if (j > 0) {
              m_Res[resIndex(i_var, i, j-1, k)] -= flux.var[i_var];
              countFlops(1);
            }
            m_Res[resIndex(i_var, i, j, k)]   += flux.var[i_var];
          }
          countFlops(DIM);
        }

        // z direction
        if (k > 0) {
          flux.fill(0);
          m_Flux.z(patch(), i, j, k, Az, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            if (k > 0) {
              m_Res[resIndex(i_var, i, j, k-1)] -= flux.var[i_var];
              countFlops(1);
            }
            m_Res[resIndex(i_var, i, j, k)]   += flux.var[i_var];
          }
          countFlops(DIM);
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
        flux.fill(0);
        m_Flux.xWallM(patch(), i1, j, k, Ax, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i1, j, k)] += flux.var[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // .. right wall
  if (i2 == patch()->sizeI()) {
    for (size_t j = j1; j < j2; ++j) {
      for (size_t k = k1; k < k2; ++k) {
        flux.fill(0);
        m_Flux.xWallP(patch(), i2-1, j, k, Ax, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i2-1, j, k)] -= flux.var[i_var];
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
        flux.fill(0);
        m_Flux.yWallM(patch(), i, j1, k, Ay, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j1, k)] += flux.var[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // .. back wall
  if (j2 == patch()->sizeJ()) {
    for (size_t i = i1; i < i2; ++i) {
      for (size_t k = k1; k < k2; ++k) {
        flux.fill(0);
        m_Flux.yWallP(patch(), i, j2-1, k, Ay, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j2-1, k)] -= flux.var[i_var];
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
        flux.fill(0);
        m_Flux.zWallM(patch(), i, j, k1, Az, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j, k1)] += flux.var[i_var];
        }
        countFlops(DIM);
      }
    }
  }

  // .. top wall
  if (k2 == patch()->sizeK()) {
    for (size_t i = i1; i < i2; ++i) {
      for (size_t j = j1; j < j2; ++j) {
        flux.fill(0);
        m_Flux.zWallP(patch(), i, j, k2-1, Az, flux);
        for (size_t i_var = 0; i_var < DIM; ++i_var) {
          m_Res[resIndex(i_var, i, j, k2-1)] -= flux.var[i_var];
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

#endif // CARTESIANSTANDARDPATCHOPERATION_H
