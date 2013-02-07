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


void StructuredHexRaster::resize(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  m_NumL = m_NumI * m_NumJ * m_NumK; // = 1
//  computeDeltas();
//  if(m_CellLink) {
//    delete m_CellLink;
//  }
//  size_t num_cells = m_NumI * m_NumJ * m_NumK;
//  m_CellLink = new List(num_cells, num_cells/10);  /// @todo Our version of "List" sets a minimum increment
}


