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

private: // methods



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
   * @brief Get the value of a field at an (i, j, k) triple.
   * @param field a pointer to the data field
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the field value at (i, j, k).
   */
  real& f(real *var, int i, int j, int k) { return var[i*m_NumJ*m_NumK + j*m_NumK + k]; }

  real rx(real *field, int i1, int i2, int j, int k); ///< gradient ratio in x direction for limiters
  real ry(real *field, int i, int j1, int j2, int k); ///< gradient ratio in y direction for limiters
  real rz(real *field, int i, int j, int k1, int k2); ///< gradient ratio in z direction for limiters

  real dx() { return m_DX; }
  real dy() { return m_DY; }
  real dz() { return m_DZ; }
  real idx() { return m_InvDX; }
  real idy() { return m_InvDY; }
  real idz() { return m_InvDZ; }

  real gasR()     { return 287; }    ///< @todo find a concept for this
  real gasGamma() { return 1.4; }    ///< @todo find a concept for this
  real gasCp()    { return 1004.5; } ///< @todo find a concept for this
  real gasCv()    { return 717.5; }  ///< @todo find a concept for this

#ifdef WITH_VTK
  void writeToVtk(QString file_name);
#endif

};

inline real CartesianPatch::rx(real *field, int i1, int i2, int j, int k)
{
  return (f(field, i1, j, k) - f(field, 2*i1 - i2, j, k))/nonZero(f(field, i2, j, k) - f(field, i1, j, k), m_LimiterEpsilon);
}

inline real CartesianPatch::ry(real *field, int i, int j1, int j2, int k)
{
  return (f(field, i, j1, k) - f(field, i, 2*j1 - j2, k))/nonZero(f(field, i, j2, k) - f(field, i, j1, k), m_LimiterEpsilon);
}

inline real CartesianPatch::rz(real *field, int i, int j, int k1, int k2)
{
  return (f(field, i, j, k1) - f(field, i, j, 2*k1 - k2))/nonZero(f(field, i, j, k2) - f(field, i, j, k1), m_LimiterEpsilon);
}


#define PROJ_UX(F, I1, I2,  J,  K) \
  f(F, I1,  J,  K) + 0.5* LIM(rx(F, I1, I2,  J,  K)) * (f(F, I1,  J,  K) - f(F, 2*I1-I2,  J,  K))

#define PROJ_UY(F,  I, J1, J2,  K) \
  f(F,  I, J1,  K) + 0.5* LIM(ry(F,  I, J1, J2,  K)) * (f(F,  I, J1,  K) - f(F,  I, 2*J1-J2,  K))

#define PROJ_UZ(F,  I,  J, K1, K2) \
  f(F,  I,  J, K1) + 0.5* LIM(ry(F,  I,  J, K1, K2)) * (f(F,  I,  J, K1) - f(F,  I,  J, 2*K1-J2))


#endif // CARTESIANPATCH_H
