#include "structuredhexraster.h"

StructuredHexRaster::StructuredHexRaster()
{
  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
  m_NumL = m_NumI * m_NumJ * m_NumK; // = 1
  //m_CellLink = new List(1,1);
  //m_CellLink = NULL;
  m_Eps = 1.e-6;  /// @todo need a better eps-handling.
}


void StructuredHexRaster::resize(const size_t& num_i, const size_t& num_j, const size_t& num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  m_NumL = m_NumI * m_NumJ * m_NumK; // = 1
}


