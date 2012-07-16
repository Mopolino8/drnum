#include "cartesianpatch.h"

CartesianPatch::CartesianPatch()
{
}

void CartesianPatch::allocateData(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  Patch::allocateData(m_NumI*m_NumJ*m_NumK);
}
