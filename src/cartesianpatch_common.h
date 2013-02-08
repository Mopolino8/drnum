#include "blockcfd_cuda.h"

private:

size_t m_NumI;
size_t m_NumJ;
size_t m_NumK;

real m_Dx;
real m_Dy;
real m_Dz;
real m_InvDx;
real m_InvDy;
real m_InvDz;


protected:


public:

CUDA_DH size_t sizeI() { return m_NumI; }
CUDA_DH size_t sizeJ() { return m_NumJ; }
CUDA_DH size_t sizeK() { return m_NumK; }

CUDA_DH real dx()  { return m_Dx; }
CUDA_DH real dy()  { return m_Dy; }
CUDA_DH real dz()  { return m_Dz; }
CUDA_DH real dV()  { return m_Dx*m_Dy*m_Dz; }
CUDA_DH real idx() { return m_InvDx; }
CUDA_DH real idy() { return m_InvDy; }
CUDA_DH real idz() { return m_InvDz; }

/**
 * @brief Get the field index of an (i, j, k) triple/
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @return the index in the one dimensional data field
 */
CUDA_DH size_t index(int i, int j, int k) { return i*m_NumJ*m_NumK + j*m_NumK + k; }

/**
 * @brief Get a variable set at a specified (i,j,k) position.
 * @param i_field the field index
 * @param i first Cartesian index
 * @param j second Cartesian index
 * @param k third Cartesian index
 * @param var will hold the conservative variable set afterwards (needs to be allocated beforehand)
 */
CUDA_DH void getVar(size_t i_field, size_t i, size_t j, size_t k, real* var)
{
  GlobalDebug::ijk(i, j, k);
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
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
CUDA_DH void setVar(size_t i_field, size_t i, size_t j, size_t k, real* var)
{
  GlobalDebug::ijk(i, j, k);
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
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
  return getVariable(i_field, i_var)[i*m_NumJ*m_NumK + j*m_NumK + k];
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
CUDA_DH void getXGrad(size_t i_field, size_t i, size_t j, size_t k, real* grad)
{
  real* var = new real[numVariables()];
  real D = 1.0/dx();
  countFlops(1);
  if (i > 0 && i < sizeI()-1) {
    getVar(i_field, i+1, j, k, grad);
    getVar(i_field, i-1, j, k, var);
    D *= 0.5;
    countFlops(1);
  } else if (i > 0) {
    getVar(i_field, i, j, k, grad);
    getVar(i_field, i-1, j, k, var);
  } else {
    getVar(i_field, i+1, j, k, grad);
    getVar(i_field, i, j, k, var);
  }
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    grad[i_var] -= var[i_var];
    grad[i_var] *= D;
  }
  countFlops(2*numVariables());
  delete [] var;
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
CUDA_DH void getYGrad(size_t i_field, size_t i, size_t j, size_t k, real* grad)
{
  if (sizeJ() > 2) {
    real* var = new real[numVariables()];
    real D = 1.0/dy();
    countFlops(1);
    if (j > 0 && j < sizeJ()-1) {
      getVar(i_field, i, j+1, k, grad);
      getVar(i_field, i, j-1, k, var);
      D *= 0.5;
      countFlops(1);
    } else if (j > 0) {
      getVar(i_field, i, j, k, grad);
      getVar(i_field, i, j-1, k, var);
    } else {
      getVar(i_field, i, j+1, k, grad);
      getVar(i_field, i, j, k, var);
    }
    for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
      grad[i_var] -= var[i_var];
      grad[i_var] *= D;
    }
    countFlops(2*numVariables());
    delete [] var;
  } else {
    for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
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
CUDA_DH void getZGrad(size_t i_field, size_t i, size_t j, size_t k, real* grad)
{
  real* var = new real[numVariables()];
  real D = 1.0/dz();
  countFlops(1);
  if (k > 0 && k < sizeK()-1) {
    getVar(i_field, i, j, k+1, grad);
    getVar(i_field, i, j, k-1, var);
    D *= 0.5;
    countFlops(1);
  } else if (k > 0) {
    getVar(i_field, i, j, k, grad);
    getVar(i_field, i, j, k-1, var);
  } else {
    getVar(i_field, i, j, k+1, grad);
    getVar(i_field, i, j, k, var);
  }
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    grad[i_var] -= var[i_var];
    grad[i_var] *= D;
  }
  countFlops(2*numVariables());
  delete [] var;
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
  GlobalDebug::ijk(i, j, k);
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
CUDA_HO void copyAttributes(const T& obj)
{
  m_NumI  = obj.sizeI();
  m_NumJ  = obj.sizeJ();
  m_NumK  = obj.sizeK();
  m_Dx    = obj.dx();
  m_Dy    = obj.dy();
  m_Dz    = obj.dz();
  m_InvDx = obj.idx();
  m_InvDy = obj.idy();
  m_InvDz = obj.idz();
}







