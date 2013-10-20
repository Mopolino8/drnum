#ifndef LSLAYERDATA_H
#define LSLAYERDATA_H

/**
  * Data sets to be stored per layer cells.
  *
  * Cell index
  * Levelset value
  * Levelset gradient in xo,yo,zo coordinates
  */
struct LSLayerData
{
  size_t m_Cell; // index of the cell (in a patch)
  real m_G;      // Levelset value
  real m_Gx;     // x-comp of gradient
  real m_Gy;     // y-comp of gradient
  real m_Gz;     // z-comp of gradient

  LSLayerData(){}

  LSLayerData(size_t cell) {
    m_Cell = cell;
  }

  LSLayerData(size_t cell, real g) {
    m_Cell = cell;
    m_G = g;
  }

  void operator= (const LSLayerData& other_ls_layer_data) {
    m_Cell = other_ls_layer_data.m_Cell;
    m_G    = other_ls_layer_data.m_G;
    m_Gx   = other_ls_layer_data.m_Gx;
    m_Gy   = other_ls_layer_data.m_Gy;
    m_Gz   = other_ls_layer_data.m_Gz;
  }

};

#endif // LSLAYERDATA_H
