#include "compressibleeulerlsobc.h"


CompressibleEulerLSOBC::CompressibleEulerLSOBC (size_t field, LevelSetObject* levelset_object)
  : LevelSetObjectBC (field, levelset_object)
{
}


void CompressibleEulerLSOBC::operator()()
{
  // ATTENTION Assume first 5 variables to make up compressible vars
  // in the sequence rho, rhou, rhov, rhow, rhoE

  // Apply mirror-T-wall condition on all cells in given inner cell lists

  real relax = 0.5;
  /** @todo a better relaxation system would help: Reduce/omit unphysical influence of
    *       relax < 1. , if corrections are small, but keep it to react on inpulsive
    *       starts or similar. Found approx. 0.8 to be a limit for impulsive starts. */

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

    llc_start = m_InnerCLStart[i_p][0];
    //TEST llc_start = m_InnerCLStart[i_p][1];


    llc_end   = m_InnerCLStart[i_p+1][0];
    for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
      /** @todo Instead of excluding action via m_ExOK, it would be better to
        *       exclude the cell from m_InnerCellsLayers. */
      if(m_InnerCellsLayers[ll_c].m_ExOK) {
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

        // Apply mirror condition;
        //  * copy rho_acc and rhoE_acc => keep p, T, ... (as also |rho_uvw| is kept)
        //  * keep tangential speed
        //  * reverse normal speed found on access point
        // rhouvw_mirror = rhouvw_acc - 2*(nabla_G * rhouvw) * nabla_G
        //
        //.. Normal projection of rhouvw_acc
        rho_uvw_n = rhou_acc * gxn + rhov_acc * gyn + rhow_acc * gzn;
        //.. Apply mirror condition
        //   Experimental version: reduce influence of relaxation, if corrections
        //   are small anyway.
        real eps = 1.e-6;
        real rhou_hard, rhov_hard, rhow_hard;
        real rhou_diff, rhov_diff, rhow_diff;
        real damp_rhou, damp_rhov, damp_rhow;
        rhou_hard = rhou_acc - 2. * rho_uvw_n * gxn;
        rhov_hard = rhov_acc - 2. * rho_uvw_n * gyn;
        rhow_hard = rhow_acc - 2. * rho_uvw_n * gzn;
        rhou_diff = rhou[l_c] - rhou_hard;
        rhov_diff = rhov[l_c] - rhov_hard;
        rhow_diff = rhow[l_c] - rhow_hard;
        damp_rhou = abs(rhou_diff) / (abs(rhou[l_c]) + eps);
        damp_rhov = abs(rhov_diff) / (abs(rhov[l_c]) + eps);
        damp_rhow = abs(rhow_diff) / (abs(rhow[l_c]) + eps);
        if(damp_rhou > (1.-relax)) damp_rhou = 1. - relax;
        if(damp_rhov > (1.-relax)) damp_rhov = 1. - relax;
        if(damp_rhow > (1.-relax)) damp_rhow = 1. - relax;
        rhou[l_c] = (1. - damp_rhou) * rhou_hard + damp_rhou * rhou[l_c];
        rhov[l_c] = (1. - damp_rhov) * rhov_hard + damp_rhov * rhov[l_c];
        rhow[l_c] = (1. - damp_rhow) * rhow_hard + damp_rhow * rhow[l_c];
        // End eperimental version
        // Classical version: constant relaxation factor
//        rhou[l_c] = relax * (rhou_acc - 2. * rho_uvw_n * gxn) + (1. -relax) * rhou[l_c];
//        rhov[l_c] = relax * (rhov_acc - 2. * rho_uvw_n * gyn) + (1. -relax) * rhov[l_c];
//        rhow[l_c] = relax * (rhow_acc - 2. * rho_uvw_n * gzn) + (1. -relax) * rhow[l_c];
        // End classical version: constant relaxation factor

        //.. other variables
        rho[l_c]  = rho_acc;
        rhoE[l_c] = rhoE_acc;

        //        // Apply T-Wall condition
        //        //  * set normal speed component to zero
        //        //s
        //        //.. copy onto inner cell
        //        rho[l_c]  = rho_acc;
        //        rhou[l_c] = rhou_acc;
        //        rhov[l_c] = rhov_acc;
        //        rhow[l_c] = rhow_acc;
        //        rhoE[l_c] = rhoE_acc;
        //        //.. apply T-Wall condition (set normal speed to 0)
        //        rho_uvw_n = rhou[l_c] * gxn + rhov[l_c] * gyn + rhow[l_c] * gzn;
        //        rhou[l_c] -= rho_uvw_n * gxn;
        //        rhov[l_c] -= rho_uvw_n * gyn;
        //        rhow[l_c] -= rho_uvw_n * gzn;
      }
    }

    // Outer cells;  Do nothing

    //    llc_start = m_OuterCLStart[i_p][0];
    //    llc_end   = m_OuterCLStart[i_p+1][0];
    //    for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
    //      l_c = m_OuterCellsLayers[ll_c].m_Cell;
    //      gxn = m_OuterCellsLayers[ll_c].m_Gx;
    //      gyn = m_OuterCellsLayers[ll_c].m_Gy;
    //      gzn = m_OuterCellsLayers[ll_c].m_Gz;
    //      rho_uvw_n = rhou[l_c] * gxn + rhov[l_c] * gyn + rhow[l_c] * gzn;
    //      rhou[l_c] -= rho_uvw_n * gxn;
    //      rhov[l_c] -= rho_uvw_n * gyn;
    //      rhow[l_c] -= rho_uvw_n * gzn;
    //    }
  }
}
