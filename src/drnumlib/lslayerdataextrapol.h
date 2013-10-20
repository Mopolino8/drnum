#ifndef LSLAYERDATAEXTRAPOL_H
#define LSLAYERDATAEXTRAPOL_H

/**
  * Data sets to be stored per layer cells.
  *
  * Data content:
  *   Cell index
  *   Levelset value
  *   Levelset gradient in xo,yo,zo coordinates
  *   8 pairs of indices and weights for extrapolation data access
  *
  * Note: At most 8 interpolation partners are allowed. In case, more
  *       partners are present (unstructured grids), the 8 most dominant
  *       must be selected.
  *
  */
struct LSLayerDataExtrapol
{
  size_t m_Cell;       /// index of the cell (in a patch)

  real m_G;            /// Levelset value
  real m_Gx;           /// x-comp of gradient
  real m_Gy;           /// y-comp of gradient
  real m_Gz;           /// z-comp of gradient

  bool m_ExOK;         /// flag, true if a valid extrapolate is available

  size_t m_MirrorDonor[8]; /// indices for mirror interpolation donors
  real   m_MirrorWeight[8];  /// weights for mirror interpolation donors


  LSLayerDataExtrapol() {
    m_ExOK = false;
  }

  LSLayerDataExtrapol(size_t cell) {
    m_Cell = cell;
    m_ExOK = false;
  }

  LSLayerDataExtrapol(size_t cell, real g) {
    m_Cell = cell;
    m_G = g;
    m_ExOK = false;
  }

  void operator= (const LSLayerDataExtrapol& other_ls_layer_data) {
    m_Cell = other_ls_layer_data.m_Cell;
    m_G    = other_ls_layer_data.m_G;
    m_Gx   = other_ls_layer_data.m_Gx;
    m_Gy   = other_ls_layer_data.m_Gy;
    m_Gz   = other_ls_layer_data.m_Gz;
    m_ExOK = other_ls_layer_data.m_ExOK;
    for(size_t i = 0; i<8; i++) {
      m_MirrorDonor[i] = other_ls_layer_data.m_MirrorDonor[i];
      m_MirrorWeight[i] = other_ls_layer_data.m_MirrorWeight[i];
    }
  }

  real accessInterpolate(real* var) {
    return (
          m_MirrorWeight[0] * var[m_MirrorDonor[0]] +
          m_MirrorWeight[1] * var[m_MirrorDonor[1]] +
          m_MirrorWeight[2] * var[m_MirrorDonor[2]] +
          m_MirrorWeight[3] * var[m_MirrorDonor[3]] +
          m_MirrorWeight[4] * var[m_MirrorDonor[4]] +
          m_MirrorWeight[5] * var[m_MirrorDonor[5]] +
          m_MirrorWeight[6] * var[m_MirrorDonor[6]] +
          m_MirrorWeight[7] * var[m_MirrorDonor[7]]
          );
  }

};

#endif // LSLAYERDATAEXTRAPOL_H
