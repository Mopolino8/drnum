#include "compressibleeulerlsobc.h"


CompressibleEulerLSOBC::CompressibleEulerLSOBC (size_t field, LevelSetObject* levelset_object)
  : LevelSetObjectBC (field, levelset_object)
{
}


void CompressibleEulerLSOBC::operator()()
{
  // ATTENTION Assume first 5 variables to make up compressible vars
  // in the sequence rho, rhou, rhov, rhow, rhoE

  // Apply T-wall condition on all cells in given lists
  size_t num_patches = m_PatchGrid->getNumPatches();
  for (size_t i_p = 0; i_p < num_patches; i_p++) {
    Patch* patch = m_LevelSetObject->getPatchGrid()->getPatch(i_p);
    real* rho  = patch->getVariable(m_Field, 0);  // 0: rho
    real* rhou = patch->getVariable(m_Field, 1);  // 1: rhou
    real* rhov = patch->getVariable(m_Field, 2);  // 2: rhov
    real* rhow = patch->getVariable(m_Field, 3);  // 3: rhow
    real* rhoE = patch->getVariable(m_Field, 4);  // 4: rhoE
    size_t llc_start, llc_end;
    size_t l_c;
    real gxn, gyn, gzn;
    real rho_uvw_n;

//    /// @todo DEBUG only
//    real* g = patch->getVariable(m_Field, 5);

    //.. Inner cells
    //   Note: Do all layers at once

    //llc_start = m_InnerCLStart[i_p][0];
    llc_start = m_InnerCLStart[i_p][1];


    llc_end   = m_InnerCLStart[i_p+1][0];
    for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
      l_c = m_InnerCellsLayers[ll_c].m_Cell;
      gxn = m_InnerCellsLayers[ll_c].m_Gx;
      gyn = m_InnerCellsLayers[ll_c].m_Gy;
      gzn = m_InnerCellsLayers[ll_c].m_Gz;
      //.. access mirror data
      real rho_acc = 0.;
      real rhou_acc = 0.;
      real rhov_acc = 0.;
      real rhow_acc = 0.;
      real rhoE_acc = 0.;
      for(size_t i = 0; i < 8; i++) {
        size_t donor = m_InnerCellsLayers[ll_c].m_MirrorDonor[i];
        real weight = m_InnerCellsLayers[ll_c].m_MirrorWeight[i];
        rho_acc  += weight * rho[donor];
        rhou_acc += weight * rhou[donor];
        rhov_acc += weight * rhov[donor];
        rhow_acc += weight * rhow[donor];
        rhoE_acc += weight * rhoE[donor];
      }
      //.. copy onto inner cell
      /// @todo needs management on access failure, see m_InnerCellsLayers[ll_c].m_ExOK
      rho[l_c]  = rho_acc;
      rhou[l_c] = rhou_acc;
      rhov[l_c] = rhov_acc;
      rhow[l_c] = rhow_acc;
      rhoE[l_c] = rhoE_acc;

      //.. apply T-Wall condition (set normal speed to 0)
      rho_uvw_n = rhou[l_c] * gxn + rhov[l_c] * gyn + rhow[l_c] * gzn;
      rhou[l_c] -= rho_uvw_n * gxn;
      rhov[l_c] -= rho_uvw_n * gyn;
      rhow[l_c] -= rho_uvw_n * gzn;
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

//      /// @todo DEBUG only
//      real x_c, y_c, z_c;
//      patch->xyzCell(l_c,
//                     x_c, y_c, z_c);
//      real gc = g[l_c];
//      gc = g[l_c];

    }
  }
}
