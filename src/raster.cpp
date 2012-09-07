#include "raster.h"

Raster::Raster()
{
  //m_CellLink = new List(1,1);
  //m_CellLink = NULL;
  //m_Eps = 1.e-6;  /// @todo need a better eps-handling.
}


void Raster::resize(size_t num_l)
{
  m_NumL = num_l;
}


void Raster::setNumProtectLayers(size_t num_protectlayers)
{
  m_NumProtectLayers = num_protectlayers;
}
