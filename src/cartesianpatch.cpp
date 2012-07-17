#include "cartesianpatch.h"

CartesianPatch::CartesianPatch()
{
  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
}

void CartesianPatch::allocateData(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  Patch::allocateData(m_NumI*m_NumJ*m_NumK);
}

void CartesianPatch::computeDeltas()
{
  double Lx = sqrt(sqr(m_UX) + sqr(m_UY) + sqr(m_UZ));
  double Ly = sqrt(sqr(m_VX) + sqr(m_VY) + sqr(m_VZ));
  double Lz = sqrt(sqr(m_WX) + sqr(m_WY) + sqr(m_WZ));

  m_DX = Lx/m_NumI;
  m_DY = Ly/m_NumJ;
  m_DZ = Lz/m_NumK;

  m_InvDX = 1.0/m_DX;
  m_InvDY = 1.0/m_DX;
  m_InvDZ = 1.0/m_DX;
}

void CartesianPatch::setupAligned(real x1, real y1, real z1, real x2, real y2, real z2)
{
  m_X0 = x1;
  m_Y0 = y1;
  m_Z0 = z1;

  m_UX = x2 - x1;
  m_UY = 0;
  m_UZ = 0;

  m_VX = 0;
  m_VY = y2 = y1;
  m_VZ = 0;

  m_WX = 0;
  m_WY = 0;
  m_WZ = z2 - z1;

  computeDeltas();
}

void CartesianPatch::resize(size_t num_i, size_t num_j, size_t num_k)
{
  /// @todo implement a field mapping mechanism if the grid gets resized
  deleteData();
  allocateData(num_i, num_j, num_k);
  computeDeltas();
}
