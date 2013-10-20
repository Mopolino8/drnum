#include "compressiblesimpleswalllsobc.h"


CompressibleSimpleSWallLSOBC::CompressibleSimpleSWallLSOBC (size_t field, LevelSetObject* levelset_object)
  : LevelSetObjectBC (field, levelset_object)
{
}


void CompressibleSimpleSWallLSOBC::operator()()
{
  // ATTENTION Assume first 5 variables to make up compressible vars
  // in the sequence rho, rhou, rhov, rhow, rhoE

  // Apply:
  //  * simple S-wall condition on all cells in given inner lists
  //  * T-wall condition on all cells in given outer lists

  size_t num_patches = m_PatchGrid->getNumPatches();
  for (size_t i_p = 0; i_p < num_patches; i_p++) {
    Patch* patch = m_LevelSetObject->getPatchGrid()->getPatch(i_p);
    real* rhou = patch->getVariable(m_Field, 1);  // 1: rhou
    real* rhov = patch->getVariable(m_Field, 2);  // 2: rhov
    real* rhow = patch->getVariable(m_Field, 3);  // 3: rhow
    size_t llc_start, llc_end;
    size_t l_c;
    real gxn, gyn, gzn;
    real rho_uvw_n;

    //.. Inner cells
    //   Note: Do all layers at once
    llc_start = m_InnerCLStart[i_p][0];
    llc_end   = m_InnerCLStart[i_p+1][0];
    for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
      l_c = m_InnerCellsLayers[ll_c].m_Cell;
      rhou[l_c] = 0.;
      rhov[l_c] = 0.;
      rhow[l_c] = 0.;
    }

    //.. Outer cells
    //   Note: Do all layers at once

    llc_start = m_OuterCLStart[i_p][0];
    llc_end   = m_OuterCLStart[i_p+1][0];
    for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
      l_c = m_OuterCellsLayers[ll_c].m_Cell;
      gxn = m_OuterCellsLayers[ll_c].m_Gx;
      gyn = m_OuterCellsLayers[ll_c].m_Gy;
      gzn = m_OuterCellsLayers[ll_c].m_Gz;
      rho_uvw_n = rhou[l_c] * gxn + rhov[l_c] * gyn + rhow[l_c] * gzn;
      rhou[l_c] -= rho_uvw_n * gxn;
      rhov[l_c] -= rho_uvw_n * gyn;
      rhow[l_c] -= rho_uvw_n * gzn;
    }

  }
}
