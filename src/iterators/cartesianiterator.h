#ifndef CARTESIANITERATOR_H
#define CARTESIANITERATOR_H

#include "cartesianpatch.h"
#include "iterators/tpatchiterator.h"

template <unsigned int DIM, typename OP>
class CartesianIterator : public TPatchIterator<CartesianPatch, DIM, OP>
{

protected: // attributes

  real*  m_Res;
  size_t m_ResLength;
  size_t m_I1;
  size_t m_J1;
  size_t m_K1;
  size_t m_SizeI;
  size_t m_SizeJ;
  size_t m_SizeK;


protected: // methods

  void checkResFieldSize(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2, size_t num_vars);


public:

  using TPatchIterator<CartesianPatch, DIM, OP>::addPatch;

  CartesianIterator(PatchGrid &patch_grid, OP op);

  size_t resIndex(size_t i_var, size_t i, size_t j, size_t k) { return m_ResLength*i_var + (i-m_I1)*m_SizeJ*m_SizeK + (j-m_J1)*m_SizeK + (k-m_K1); }

  virtual void compute(real factor, const vector<size_t> &patches);

};

template <unsigned int DIM, typename OP>
CartesianIterator<DIM, OP>::CartesianIterator(PatchGrid &patch_grid, OP op) : TPatchIterator<CartesianPatch, DIM, OP>(patch_grid, op)
{
  m_Res = NULL;
  m_ResLength = 0;
}

template <unsigned int DIM, typename OP>
void CartesianIterator<DIM, OP>::checkResFieldSize(size_t i1, size_t j1, size_t k1, size_t i2, size_t j2, size_t k2, size_t num_vars)
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
    m_Res = new real [new_length*num_vars];
    m_ResLength = new_length;
  }
}


template <unsigned int DIM, typename OP>
void CartesianIterator<DIM, OP>::compute(real factor, const vector<size_t> &patches)
{
  for (size_t i_patch = 0; i_patch < patches.size(); ++i_patch) {
    CartesianPatch* patch = this->m_Patches[patches[i_patch]];

    checkResFieldSize(0, 0, 0, patch->sizeI(), patch->sizeJ(), patch->sizeK(), patch->numVariables());

    real Ax = patch->dy()*patch->dz();
    real Ay = patch->dx()*patch->dz();
    real Az = patch->dx()*patch->dy();

    size_t i1 = 0;
    size_t j1 = 0;
    size_t k1 = 0;
    size_t i2 = patch->sizeI();
    size_t j2 = patch->sizeJ();
    size_t k2 = patch->sizeK();

    real flux[5];

    for (size_t i_res = 0; i_res < DIM*m_ResLength; ++i_res) {
      m_Res[i_res] = 0;
    }
    countFlops(3);

    // compute main block
    for (int offset = 0; offset <= 1; ++offset) {
      #ifndef DEBUG
      #pragma omp parallel
      #endif
      {
        size_t num_threads = omp_get_num_threads();
        size_t tid         = omp_get_thread_num();
        size_t n           = patch->sizeI()/(2*num_threads) + 1;
        size_t i_start     = max(i1, (offset + 2*tid)*n);
        size_t i_stop      = min(i2, i_start + n);

        real flux[5];

        #ifdef DEBUG
        if (num_threads != 1) {
          BUG;
        }
        #endif

        real x = 0.5*patch->dx() + i_start * patch->dx();
        for (size_t i = i_start; i < i_stop; ++i) {
          real y = 0.5*patch->dy();
          for (size_t j = j1; j < j2; ++j) {
            real z = 0.5*patch->dz();
            for (size_t k = k1; k < k2; ++k) {

              GlobalDebug::xyz(x,y,z);

              // x direction
              if (i > 0 && patch->sizeI() > 2) {
                fill(flux, 5, 0);
                this->m_Op.xField(patch, i, j, k, x, y, z, Ax, flux);
                for (size_t i_var = 0; i_var < DIM; ++i_var) {
                  m_Res[resIndex(i_var, i-1, j, k)] -= flux[i_var];
                  m_Res[resIndex(i_var, i, j, k)]   += flux[i_var];
                }
                countFlops(2*DIM);
              }

              // y direction
              if (j > 0 && patch->sizeJ() > 2) {
                fill(flux, 5, 0);
                this->m_Op.yField(patch, i, j, k, x, y, z, Ay, flux);
                for (size_t i_var = 0; i_var < DIM; ++i_var) {
                  m_Res[resIndex(i_var, i, j-1, k)] -= flux[i_var];
                  m_Res[resIndex(i_var, i, j, k)]   += flux[i_var];
                }
                countFlops(2*DIM);
              }

              // z direction
              if (k > 0 && patch->sizeK() > 2) {
                fill(flux, 5, 0);
                this->m_Op.zField(patch, i, j, k, x, y, z, Az, flux);
                for (size_t i_var = 0; i_var < DIM; ++i_var) {
                  m_Res[resIndex(i_var, i, j, k-1)] -= flux[i_var];
                  m_Res[resIndex(i_var, i, j, k)]   += flux[i_var];
                }
                countFlops(2*DIM);
              }

              z += patch->dz();
            }
            y += patch->dy();
          }
          x += patch->dx();
        }
      }
    }

    // compute x walls
    //
    // .. left wall
    if (i1 == 0 && patch->sizeI() > 2) {
      real x = 0.5*patch->dx();
      real y = 0.5*patch->dy();
      for (size_t j = j1; j < j2; ++j) {
        real z = 0.5*patch->dz();
        for (size_t k = k1; k < k2; ++k) {
          fill(flux, 5, 0);
          this->m_Op.xWallM(patch, i1, j, k, x, y, z, Ax, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i1, j, k)] += flux[i_var];
          }
          countFlops(DIM);
          z += patch->dz();
        }
        y += patch->dy();
      }
    }

    // .. right wall
    if (i2 == patch->sizeI() && patch->sizeI() > 2) {
      real x = 0.5*patch->dx() + patch->sizeI()*patch->dx();
      real y = 0.5*patch->dy();
      for (size_t j = j1; j < j2; ++j) {
        real z = 0.5*patch->dz();
        for (size_t k = k1; k < k2; ++k) {
          fill(flux, 5, 0);
          this->m_Op.xWallP(patch, i2, j, k, x, y, z, Ax, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i2-1, j, k)] -= flux[i_var];
          }
          countFlops(DIM);
          z += patch->dz();
        }
        y += patch->dy();
      }
    }

    // compute y walls
    //
    // .. front wall
    if (j1 == 0 && patch->sizeJ() > 2) {
      real x = 0.5*patch->dx();
      real y = 0.5*patch->dy();
      for (size_t i = i1; i < i2; ++i) {
        real z = 0.5*patch->dz();
        for (size_t k = k1; k < k2; ++k) {
          fill(flux, 5, 0);
          this->m_Op.yWallM(patch, i, j1, k, x, y, z, Ay, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i, j1, k)] += flux[i_var];
          }
          countFlops(DIM);
          z += patch->dz();
        }
        x += patch->dx();
      }
    }

    // .. back wall
    if (j2 == patch->sizeJ() && patch->sizeJ() > 2) {
      real x = 0.5*patch->dx();
      real y = 0.5*patch->dy() + patch->sizeJ()*patch->dy();
      for (size_t i = i1; i < i2; ++i) {
        real z = 0.5*patch->dz();
        for (size_t k = k1; k < k2; ++k) {
          fill(flux, 5, 0);
          this->m_Op.yWallP(patch, i, j2, k, x, y, z, Ay, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i, j2-1, k)] -= flux[i_var];
          }
          countFlops(DIM);
          z += patch->dz();
        }
        x += patch->dx();
      }
    }

    // compute z walls
    //
    // .. bottom wall
    if (k1 == 0 && patch->sizeK() > 2) {
      real x = 0.5*patch->dx();
      real z = 0.5*patch->dz();
      for (size_t i = i1; i < i2; ++i) {
        real y = 0.5*patch->dy();
        for (size_t j = j1; j < j2; ++j) {
          fill(flux, 5, 0);
          this->m_Op.zWallM(patch, i, j, k1, x, y, z, Az, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i, j, k1)] += flux[i_var];
          }
          countFlops(DIM);
          y += patch->dy();
        }
        x += patch->dx();
      }
    }

    // .. top wall
    if (k2 == patch->sizeK() && patch->sizeK() > 2) {
      real x = 0.5*patch->dx();
      real z = 0.5*patch->dz() + patch->sizeK()*patch->dz();
      for (size_t i = i1; i < i2; ++i) {
        real y = 0.5*patch->dy();
        for (size_t j = j1; j < j2; ++j) {
          fill(flux, 5, 0);
          this->m_Op.zWallP(patch, i, j, k2, x, y, z, Az, flux);
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            m_Res[resIndex(i_var, i, j, k2-1)] -= flux[i_var];
          }
          countFlops(DIM);
          y += patch->dy();
        }
        x += patch->dx();
      }
    }

    // advance to next iteration level (time)
    factor /= patch->dV();
    for (size_t i = i1; i < i2; ++i) {
      for (size_t j = j1; j < j2; ++j) {
        for (size_t k = k1; k < k2; ++k) {
          for (size_t i_var = 0; i_var < DIM; ++i_var) {
            patch->f(0, i_var, i, j, k) = patch->f(1, i_var, i, j, k) + factor*m_Res[resIndex(i_var, i, j, k)];
          }
        }
      }
    }
  }
}

#endif // CARTESIANITERATOR_H
