#include "compressibleeulerlsobc.h"


CompressibleEulerLSOBC::CompressibleEulerLSOBC (size_t field,
                                                LevelSetObject* levelset_object,
                                                size_t abuse_field)
  : LevelSetObjectBC (field, levelset_object, abuse_field)
{
}


void CompressibleEulerLSOBC::operator()()
{
  // ATTENTION Assume first 5 variables to make up compressible vars
  // in the sequence rho, rhou, rhov, rhow, rhoE

  // Apply mirror-T-wall condition on all cells in given inner cell lists

  // Potential recursion: Interpolate sets may contain cells in m_InnerCellsLayers.
  // To prevent recursion, the abuse_field is used as intermediate data storage.

  real relax = 0.5;
  /** @todo a better relaxation method would help: Reduce/omit unphysical influence of
    *       relax < 1. , if corrections are small, but keep it to react on inpulsive
    *       starts or similar. Found approx. 0.8 to be a limit for impulsive starts. */

  size_t llc_start, llc_end;
  size_t l_c;
  real gxn, gyn, gzn;
  real rho_uvw_n;

  size_t num_patches = m_PatchGrid->getNumPatches();
  llc_start = m_InnerCLStartAll[0];         // first patch
  llc_end = m_InnerCLStartAll[num_patches]; // after last patch

  // 1st loop: acquire data, avoid recursion
  for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
    /** @todo Instead of excluding action via m_ExOK, it would be better to
      *       exclude the cell from m_InnerCellsLayers. */
    if(m_InnerCellsLayers[ll_c].m_ExOK) {
      LSLayerDataExtrapol lslde_h = m_InnerCellsLayers[ll_c];
      l_c = lslde_h.m_Cell;
      //.. pointers
      real* rho   = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 0 * lslde_h.m_VariableSize;  // 0: rho
      real* rhou  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 1 * lslde_h.m_VariableSize;  // 1: rhou
      real* rhov  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 2 * lslde_h.m_VariableSize;  // 1: rhov
      real* rhow  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 3 * lslde_h.m_VariableSize;  // 1: rhow
      real* rhoE  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 4 * lslde_h.m_VariableSize;  // 1: rhoE
      real* rho_acc   = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 0 * lslde_h.m_VariableSize;
      real* rhou_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 1 * lslde_h.m_VariableSize;
      real* rhov_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 2 * lslde_h.m_VariableSize;
      real* rhow_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 3 * lslde_h.m_VariableSize;
      real* rhoE_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 4 * lslde_h.m_VariableSize;
      //      real* rho  = patch->getVariable(m_Field, 0);  // 0: rho
      //      real* rhou = patch->getVariable(m_Field, 1);  // 1: rhou
      //      real* rhov = patch->getVariable(m_Field, 2);  // 2: rhov
      //      real* rhow = patch->getVariable(m_Field, 3);  // 3: rhow
      //      real* rhoE = patch->getVariable(m_Field, 4);  // 4: rhoE
      //      real* rho_acc  = patch->getVariable(m_AbuseField, 0);  // 0: rho
      //      real* rhou_acc = patch->getVariable(m_AbuseField, 1);  // 1: rhou
      //      real* rhov_acc = patch->getVariable(m_AbuseField, 2);  // 2: rhov
      //      real* rhow_acc = patch->getVariable(m_AbuseField, 3);  // 3: rhow
      //      real* rhoE_acc = patch->getVariable(m_AbuseField, 4);  // 4: rhoE
      //.. access mirror data
      rho_acc[l_c] = 0.;
      rhou_acc[l_c] = 0.;
      rhov_acc[l_c] = 0.;
      rhow_acc[l_c] = 0.;
      rhoE_acc[l_c] = 0.;
      for(size_t i = 0; i < 8; i++) {
        size_t donor = m_InnerCellsLayers[ll_c].m_MirrorDonor[i];
        real weight = m_InnerCellsLayers[ll_c].m_MirrorWeight[i];
        rho_acc[l_c]  += weight * rho[donor];
        rhou_acc[l_c] += weight * rhou[donor];
        rhov_acc[l_c] += weight * rhov[donor];
        rhow_acc[l_c] += weight * rhow[donor];
        rhoE_acc[l_c] += weight * rhoE[donor];
      }
    }
  }

  // 2nd loop: set values
  for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
    if(m_InnerCellsLayers[ll_c].m_ExOK) {
      LSLayerDataExtrapol lslde_h = m_InnerCellsLayers[ll_c];
      l_c = lslde_h.m_Cell;
      gxn = lslde_h.m_Gx;
      gyn = lslde_h.m_Gy;
      gzn = lslde_h.m_Gz;
      //.. pointers
      real* rho   = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 0 * lslde_h.m_VariableSize;  // 0: rho
      real* rhou  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 1 * lslde_h.m_VariableSize;  // 1: rhou
      real* rhov  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 2 * lslde_h.m_VariableSize;  // 1: rhov
      real* rhow  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 3 * lslde_h.m_VariableSize;  // 1: rhow
      real* rhoE  = lslde_h.m_Data + m_Field * lslde_h.m_FieldSize + 4 * lslde_h.m_VariableSize;  // 1: rhoE
      real* rho_acc   = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 0 * lslde_h.m_VariableSize;
      real* rhou_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 1 * lslde_h.m_VariableSize;
      real* rhov_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 2 * lslde_h.m_VariableSize;
      real* rhow_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 3 * lslde_h.m_VariableSize;
      real* rhoE_acc  = lslde_h.m_Data + m_AbuseField * lslde_h.m_FieldSize + 4 * lslde_h.m_VariableSize;

      // Apply mirror condition;
      //  * copy rho_acc and rhoE_acc => keep p, T, ... (as also |rho_uvw| is kept)
      //  * keep tangential speed
      //  * reverse normal speed found on access point
      // rhouvw_mirror = rhouvw_acc - 2*(nabla_G * rhouvw) * nabla_G
      //
      //.. Normal projection of rhouvw_acc
      rho_uvw_n = rhou_acc[l_c] * gxn + rhov_acc[l_c] * gyn + rhow_acc[l_c] * gzn;
      //.. Apply mirror condition
      //   Experimental version: reduce influence of relaxation, if corrections
      //   are small anyway.
      real eps = 1.e-6;
      real rhou_hard, rhov_hard, rhow_hard;
      real rhou_diff, rhov_diff, rhow_diff;
      real damp_rhou, damp_rhov, damp_rhow;
      rhou_hard = rhou_acc[l_c]  - 2. * rho_uvw_n * gxn;
      rhov_hard = rhov_acc[l_c]  - 2. * rho_uvw_n * gyn;
      rhow_hard = rhow_acc[l_c]  - 2. * rho_uvw_n * gzn;
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
      //        rhou[l_c] = relax * (rhou_acc[l_c]  - 2. * rho_uvw_n * gxn) + (1. -relax) * rhou[l_c];
      //        rhov[l_c] = relax * (rhov_acc[l_c]  - 2. * rho_uvw_n * gyn) + (1. -relax) * rhov[l_c];
      //        rhow[l_c] = relax * (rhow_acc[l_c]  - 2. * rho_uvw_n * gzn) + (1. -relax) * rhow[l_c];
      // End classical version: constant relaxation factor

      //.. other variables
      rho[l_c]  = rho_acc[l_c] ;
      rhoE[l_c] = rhoE_acc[l_c] ;

      //        // Apply T-Wall condition
      //        //  * set normal speed component to zero
      //        //s
      //        //.. copy onto inner cell
      //        rho[l_c]  = rho_acc[l_c] ;
      //        rhou[l_c] = rhou_acc[l_c] ;
      //        rhov[l_c] = rhov_acc[l_c] ;
      //        rhow[l_c] = rhow_acc[l_c] ;
      //        rhoE[l_c] = rhoE_acc[l_c] ;
      //        //.. apply T-Wall condition (set normal speed to 0)
      //        rho_uvw_n = rhou[l_c] * gxn + rhov[l_c] * gyn + rhow[l_c] * gzn;
      //        rhou[l_c] -= rho_uvw_n * gxn;
      //        rhov[l_c] -= rho_uvw_n * gyn;
      //        rhow[l_c] -= rho_uvw_n * gzn;
    }
  }

  //  size_t num_patches = m_PatchGrid->getNumPatches();
  //  for (size_t i_p = 0; i_p < num_patches; i_p++) {
  //    Patch* patch = m_LevelSetObject->getPatchGrid()->getPatch(i_p);
  //    real* rho  = patch->getVariable(m_Field, 0);  // 0: rho
  //    real* rhou = patch->getVariable(m_Field, 1);  // 1: rhou
  //    real* rhov = patch->getVariable(m_Field, 2);  // 2: rhov
  //    real* rhow = patch->getVariable(m_Field, 3);  // 3: rhow
  //    real* rhoE = patch->getVariable(m_Field, 4);  // 4: rhoE

  //    real* rho_acc  = patch->getVariable(m_AbuseField, 0);  // 0: rho
  //    real* rhou_acc = patch->getVariable(m_AbuseField, 1);  // 1: rhou
  //    real* rhov_acc = patch->getVariable(m_AbuseField, 2);  // 2: rhov
  //    real* rhow_acc = patch->getVariable(m_AbuseField, 3);  // 3: rhow
  //    real* rhoE_acc = patch->getVariable(m_AbuseField, 4);  // 4: rhoE

  //    size_t llc_start, llc_end;
  //    size_t l_c;
  //    real gxn, gyn, gzn;
  //    real rho_uvw_n;

  //    //.. Inner cells
  //    //   Note: Do all layers at once
  //    //    llc_start = m_InnerCLStart[i_p][0];
  //    //    llc_end   = m_InnerCLStart[i_p+1][0];
  //    llc_start = m_InnerCLStartAll[i_p];
  //    llc_end = m_InnerCLStartAll[i_p+1];
  //    //.. 1st loop: acquire data, avoid recursion
  //    for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
  //      /** @todo Instead of excluding action via m_ExOK, it would be better to
  //        *       exclude the cell from m_InnerCellsLayers. */
  //      if(m_InnerCellsLayers[ll_c].m_ExOK) {
  //        l_c = m_InnerCellsLayers[ll_c].m_Cell;
  //        //.. access mirror data
  //        rho_acc[l_c] = 0.;
  //        rhou_acc[l_c] = 0.;
  //        rhov_acc[l_c] = 0.;
  //        rhow_acc[l_c] = 0.;
  //        rhoE_acc[l_c] = 0.;
  //        for(size_t i = 0; i < 8; i++) {
  //          size_t donor = m_InnerCellsLayers[ll_c].m_MirrorDonor[i];
  //          real weight = m_InnerCellsLayers[ll_c].m_MirrorWeight[i];
  //          rho_acc[l_c]  += weight * rho[donor];
  //          rhou_acc[l_c] += weight * rhou[donor];
  //          rhov_acc[l_c] += weight * rhov[donor];
  //          rhow_acc[l_c] += weight * rhow[donor];
  //          rhoE_acc[l_c] += weight * rhoE[donor];
  //        }
  //      }
  //    }
  //    //.. 2nd loop: set values
  //    for (size_t ll_c = llc_start; ll_c < llc_end; ll_c++) {
  //      if(m_InnerCellsLayers[ll_c].m_ExOK) {
  //        l_c = m_InnerCellsLayers[ll_c].m_Cell;
  //        gxn = m_InnerCellsLayers[ll_c].m_Gx;
  //        gyn = m_InnerCellsLayers[ll_c].m_Gy;
  //        gzn = m_InnerCellsLayers[ll_c].m_Gz;
  //        // Apply mirror condition;
  //        //  * copy rho_acc and rhoE_acc => keep p, T, ... (as also |rho_uvw| is kept)
  //        //  * keep tangential speed
  //        //  * reverse normal speed found on access point
  //        // rhouvw_mirror = rhouvw_acc - 2*(nabla_G * rhouvw) * nabla_G
  //        //
  //        //.. Normal projection of rhouvw_acc
  //        rho_uvw_n = rhou_acc[l_c] * gxn + rhov_acc[l_c] * gyn + rhow_acc[l_c] * gzn;
  //        //.. Apply mirror condition
  //        //   Experimental version: reduce influence of relaxation, if corrections
  //        //   are small anyway.
  //        real eps = 1.e-6;
  //        real rhou_hard, rhov_hard, rhow_hard;
  //        real rhou_diff, rhov_diff, rhow_diff;
  //        real damp_rhou, damp_rhov, damp_rhow;
  //        rhou_hard = rhou_acc[l_c]  - 2. * rho_uvw_n * gxn;
  //        rhov_hard = rhov_acc[l_c]  - 2. * rho_uvw_n * gyn;
  //        rhow_hard = rhow_acc[l_c]  - 2. * rho_uvw_n * gzn;
  //        rhou_diff = rhou[l_c] - rhou_hard;
  //        rhov_diff = rhov[l_c] - rhov_hard;
  //        rhow_diff = rhow[l_c] - rhow_hard;
  //        damp_rhou = abs(rhou_diff) / (abs(rhou[l_c]) + eps);
  //        damp_rhov = abs(rhov_diff) / (abs(rhov[l_c]) + eps);
  //        damp_rhow = abs(rhow_diff) / (abs(rhow[l_c]) + eps);
  //        if(damp_rhou > (1.-relax)) damp_rhou = 1. - relax;
  //        if(damp_rhov > (1.-relax)) damp_rhov = 1. - relax;
  //        if(damp_rhow > (1.-relax)) damp_rhow = 1. - relax;
  //        rhou[l_c] = (1. - damp_rhou) * rhou_hard + damp_rhou * rhou[l_c];
  //        rhov[l_c] = (1. - damp_rhov) * rhov_hard + damp_rhov * rhov[l_c];
  //        rhow[l_c] = (1. - damp_rhow) * rhow_hard + damp_rhow * rhow[l_c];
  //        // End eperimental version
  //        // Classical version: constant relaxation factor
  //        //        rhou[l_c] = relax * (rhou_acc[l_c]  - 2. * rho_uvw_n * gxn) + (1. -relax) * rhou[l_c];
  //        //        rhov[l_c] = relax * (rhov_acc[l_c]  - 2. * rho_uvw_n * gyn) + (1. -relax) * rhov[l_c];
  //        //        rhow[l_c] = relax * (rhow_acc[l_c]  - 2. * rho_uvw_n * gzn) + (1. -relax) * rhow[l_c];
  //        // End classical version: constant relaxation factor

  //        //.. other variables
  //        rho[l_c]  = rho_acc[l_c] ;
  //        rhoE[l_c] = rhoE_acc[l_c] ;

  //        //        // Apply T-Wall condition
  //        //        //  * set normal speed component to zero
  //        //        //s
  //        //        //.. copy onto inner cell
  //        //        rho[l_c]  = rho_acc[l_c] ;
  //        //        rhou[l_c] = rhou_acc[l_c] ;
  //        //        rhov[l_c] = rhov_acc[l_c] ;
  //        //        rhow[l_c] = rhow_acc[l_c] ;
  //        //        rhoE[l_c] = rhoE_acc[l_c] ;
  //        //        //.. apply T-Wall condition (set normal speed to 0)
  //        //        rho_uvw_n = rhou[l_c] * gxn + rhov[l_c] * gyn + rhow[l_c] * gzn;
  //        //        rhou[l_c] -= rho_uvw_n * gxn;
  //        //        rhov[l_c] -= rho_uvw_n * gyn;
  //        //        rhow[l_c] -= rho_uvw_n * gzn;
  //      }
  //    }

  // Outer cells;  Do nothing

  //    llc_start = m_OuterCLStartAll[i_p];
  //    llc_end   = m_OuterCLStartAll[i_p+1];
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
