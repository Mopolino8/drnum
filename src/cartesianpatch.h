#ifndef CARTESIANPATCH_H
#define CARTESIANPATCH_H

class CartesianPatch;

#include "patch.h"

class CartesianPatch : public Patch
{

protected: // attributes

  size_t m_NumI;
  size_t m_NumJ;
  size_t m_NumK;

  real m_X0;
  real m_Y0;
  real m_Z0;
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

protected: // methods

  void allocateData(size_t num_i, size_t num_j, size_t num_k);

public: // methods

  CartesianPatch();
  void setupAligned(real x1, real y1, real z1, real x2, real y2, real z2);
  void resize(size_t num_i, size_t num_j, size_t num_k);

};

#endif // CARTESIANPATCH_H
