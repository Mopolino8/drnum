// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "drnum.h"

private:

size_t m_NumI;
size_t m_NumJ;
size_t m_NumK;
size_t m_NumJK;

real m_Dx;
real m_Dy;
real m_Dz;
real m_InvDx;
real m_InvDy;
real m_InvDz;


public:

CUDA_DH size_t sizeI() const { return m_NumI; }
CUDA_DH size_t sizeJ() const { return m_NumJ; }
CUDA_DH size_t sizeK() const { return m_NumK; }
CUDA_DH size_t sizeL() const { return m_NumI * m_NumJ * m_NumK; }

CUDA_DH real dx()  const { return m_Dx; }
CUDA_DH real dy()  const { return m_Dy; }
CUDA_DH real dz()  const { return m_Dz; }
CUDA_DH real dV()  const { return m_Dx*m_Dy*m_Dz; }
CUDA_DH real idx() const { return m_InvDx; }
CUDA_DH real idy() const { return m_InvDy; }
CUDA_DH real idz() const { return m_InvDz; }

/**
 * @brief Get the field index of an (i, j, k) triple/
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @return the index in the one dimensional data field
 */
CUDA_DH size_t index(int i, int j, int k) { return i*m_NumJ*m_NumK + j*m_NumK + k; }

/**
 * @brief Get the node index of an (i, j, k) triple.
 * Such an index can be used as a unique identifier for a node in a Cartesian patch.
 * This index does not correspond to a "real" index within some data field!
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @return a ondimensional index
 */
CUDA_DH size_t nodeIndex(int i, int j, int k) { return i*(m_NumJ + 1)*(m_NumK + 1) + j*(m_NumK + 1) + k; }

/**
 * @brief Get the node index of an (i, j, k) triple.
 * Such an index can be used as a unique identifier for a node in a Cartesian patch.
 * This index does not correspond to a "real" index within some data field!
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @return a ondimensional index
 */
CUDA_DH size_t nodeIndex(size3_t ijk) { return ijk.i*(m_NumJ + 1)*(m_NumK + 1) + ijk.j*(m_NumK + 1) + ijk.k; }

/**
 * @brief Get a variable set at a specified (i,j,k) position.
 * @param i_field the field index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @param var will hold the conservative variable set afterwards (needs to be allocated beforehand)
 */
template <typename DIM>
CUDA_DH void getVar(DIM, size_t i_field, size_t i, size_t j, size_t k, real* var)
{
  GlobalDebug::ijk(i, j, k);
  for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
    var[i_var] = f(i_field, i_var, i, j, k);
  }
}

/**
 * @brief Get a variable set at a specified (i,j,k) position. This method take the dimension (number of variables)
 *        as real parameter and should not be used in performance critical routines.
 * @param i_field the field index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @param var will hold the conservative variable set afterwards (needs to be allocated beforehand)
 */
CUDA_DH void getVarDim(unsigned int dim, size_t i_field, size_t i, size_t j, size_t k, real* var)
{
  GlobalDebug::ijk(i, j, k);
  for (size_t i_var = 0; i_var < dim; ++i_var) {
    var[i_var] = f(i_field, i_var, i, j, k);
  }
}

/**
 * @brief Set a variable set at a specified (i,j,k) position.
 * @param i_field the field index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @param the conservative variable set
 */
template <typename DIM>
CUDA_DH void setVar(DIM, size_t i_field, size_t i, size_t j, size_t k, real* var)
{
  GlobalDebug::ijk(i, j, k);
  for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
    f(i_field, i_var, i, j, k) = var[i_var];
  }
}

/**
 * @brief Get the value of a variable at an (i, j, k) triple.
 * @param i_field field index
 * @param i_var variable index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @return the field value at (i, j, k).
 */
CUDA_DH real& f(size_t i_field, size_t i_var, size_t i, size_t j, size_t k)
{
#ifdef DEBUG
  if (!checkRange(i, j, k)) {
    BUG;
  }
#endif
  GlobalDebug::ijk(i, j, k);
  return getVariable(i_field, i_var)[i*m_NumJK + j*m_NumK + k];
}

/**
 * Get the gradient in x direction at a specifed (i,j,k) position.
 * This method will automaitically respect domain borders (one sided gradient)
 * @param i_field the field index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @param grad will hold the conservative x gradient set afterwards (needs to be allocated beforehand)
 */
template <typename DIM>
CUDA_DH void getXGrad(DIM, size_t i_field, size_t i, size_t j, size_t k, real* grad)
{
  DIM dim;

  real var[DIM::dim];
  real D = (real)0.5/dx();
  countFlops(1);
  size_t i1 = i - 1;
  size_t i2 = i + 1;
  if (i1 >= sizeI()) {
    i1 = i;
    D *= 2;
    countFlops(1);
  }
  if (i2 >= sizeI()) {
    i2 = i;
    D *= 2;
    countFlops(1);
  }
  getVar(dim, i_field, i2, j, k, grad);
  getVar(dim, i_field, i1, j, k, var);
  for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
    grad[i_var] -= var[i_var];
    grad[i_var] *= D;
  }
  countFlops(2*DIM::dim);
}

/**
 * Get the gradient in y direction at a specifed (i,j,k) position.
 * This method will automaitically respect domain borders (one sided gradient)
 * @param i_field the field index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @param grad will hold the conservative y gradient set afterwards (needs to be allocated beforehand)
 */
template <typename DIM>
CUDA_DH void getYGrad(DIM, size_t i_field, size_t i, size_t j, size_t k, real* grad)
{
  DIM dim;

  if (sizeJ() > 2) {
    real var[DIM::dim];
    real D = (real)0.5/dy();
    countFlops(1);
    size_t j1 = j - 1;
    size_t j2 = j + 1;
    if (j1 >= sizeJ()) {
      j1 = j;
      D *= 2;
      countFlops(1);
    }
    if (j2 >= sizeJ()) {
      j2 = j;
      D *= 2;
      countFlops(1);
    }
    getVar(dim, i_field, i, j2, k, grad);
    getVar(dim, i_field, i, j1, k, var);
    for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
      grad[i_var] -= var[i_var];
      grad[i_var] *= D;
    }
    countFlops(2*DIM::dim);
  } else {
    for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
      grad[i_var] = 0;
    }
  }
}

/**
 * Get the gradient in z direction at a specifed (i,j,k) position.
 * This method will automaitically respect domain borders (one sided gradient)
 * @param i_field the field index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @param grad will hold the conservative z gradient set afterwards (needs to be allocated beforehand)
 */
template <typename DIM>
CUDA_DH void getZGrad(DIM, size_t i_field, size_t i, size_t j, size_t k, real* grad)
{
  DIM dim;

  real var[DIM::dim];
  real D = (real)0.5/dz();
  countFlops(1);
  size_t k1 = k - 1;
  size_t k2 = k + 1;
  if (k1 >= sizeK()) {
    k1 = k;
    D *= 2;
    countFlops(1);
  }
  if (k2 >= sizeK()) {
    k2 = k;
    D *= 2;
    countFlops(1);
  }
  getVar(dim, i_field, i, j, k2, grad);
  getVar(dim, i_field, i, j, k1, var);
  for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
    grad[i_var] -= var[i_var];
    grad[i_var] *= D;
  }
  countFlops(2*DIM::dim);
}

/**
 * Check if an (i,j,k) triple is inside the patch.
 * Attention only the upper limit will be checked (unsigned data type).
 * @param i first index
 * @param j second index
 * @param k third index
 * @return true if it is a valid (i,j,k) triple, false otherwise
 */
CUDA_DH bool checkRange(size_t i, size_t j, size_t k)
{
  //GlobalDebug::ijk(i, j, k);
  if (i >= sizeI() || j >= sizeJ() || k >= sizeK()) {
    return false;
  }
  return true;
}

/**
 * Copy simple data attributes from another object.
 * The other object can have a different type as long as the required attributes are present.
 * param obj a constant reference to the other object
 */
template <typename T>
CUDA_HO void copyAttributes(T* obj)
{
  m_NumI  = obj->sizeI();
  m_NumJ  = obj->sizeJ();
  m_NumK  = obj->sizeK();
  m_NumJK = obj->sizeJ()*obj->sizeK();
  m_Dx    = obj->dx();
  m_Dy    = obj->dy();
  m_Dz    = obj->dz();
  m_InvDx = obj->idx();
  m_InvDy = obj->idy();
  m_InvDz = obj->idz();
}

#ifdef UNDEF_DIM
#undef DIM
#undef UNDEF_DIM
#endif


