#ifndef CARTESIANPATCH_H
#define CARTESIANPATCH_H

class CartesianPatch;

#include "patch.h"

#ifdef WITH_VTK
#include <QString>
#endif


class CartesianPatch : public Patch
{

protected: // attributes

  size_t m_NumI;
  size_t m_NumJ;
  size_t m_NumK;

  real m_Xo;
  real m_Yo;
  real m_Zo;
  real m_UX;
  real m_UY;
  real m_UZ;
  real m_VX;
  real m_VY;
  real m_VZ;
  real m_WX;
  real m_WY;
  real m_WZ;
  real m_DX;
  real m_DY;
  real m_DZ;
  real m_InvDX;
  real m_InvDY;
  real m_InvDZ;

  real m_LimiterEpsilon;

protected: // methods

  void computeDeltas();

public: // methods

  CartesianPatch();
  void setupAligned(real x1, real y1, real z1, real x2, real y2, real z2);
  void resize(size_t num_i, size_t num_j, size_t num_k);

  size_t sizeI() { return m_NumI; }
  size_t sizeJ() { return m_NumJ; }
  size_t sizeK() { return m_NumK; }

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the index in the one dimensional data field
   */
  size_t index(int i, int j, int k) { return i*m_NumJ*m_NumK + j*m_NumK + k; }

  /**
   * @brief Get the value of a variable at an (i, j, k) triple.
   * @param field a pointer to the variable data
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the field value at (i, j, k).
   */
  real& f(real *var, size_t i, size_t j, size_t k) { return var[i*m_NumJ*m_NumK + j*m_NumK + k]; }

  /**
   * @brief Get a variable set at a specified (i,j,k) position.
   * @param i_field the field index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param var will hold the conservative variable set afterwards (needs to be allocated beforehand)
   */
  void getVar(size_t i_field, size_t i, size_t j, size_t k, real* var);

  /**
   * @brief Set a variable set at a specified (i,j,k) position.
   * @param i_field the field index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param the conservative variable set
   */
  void setVar(size_t i_field, size_t i, size_t j, size_t k, real* var);

  /**
   * @brief Get the value of a variable at an (i, j, k) triple.
   * @param i_field field index
   * @param i_var variable index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the field value at (i, j, k).
   */
  real& f(size_t i_field, size_t i_var, size_t i, size_t j, size_t k) { return getVariable(i_field, i_var)[i*m_NumJ*m_NumK + j*m_NumK + k]; }

  real dx() { return m_DX; }
  real dy() { return m_DY; }
  real dz() { return m_DZ; }
  real dV() { return m_DX*m_DY*m_DZ; }
  real idx() { return m_InvDX; }
  real idy() { return m_InvDY; }
  real idz() { return m_InvDZ; }

#ifdef WITH_VTK
  void writeToVtk(QString file_name);
#endif

};


inline void CartesianPatch::getVar(size_t i_field, size_t i, size_t j, size_t k, real *var)
{
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    var[i_var] = f(i_field, i_var, i, j, k);
  }
}

inline void CartesianPatch::setVar(size_t i_field, size_t i, size_t j, size_t k, real *var)
{
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    f(i_field, i_var, i, j, k) = var[i_var];
  }
}


#endif // CARTESIANPATCH_H
