// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef LSOBCCOMPRESSIBLEEULEROP_H
#define LSOBCCOMPRESSIBLEEULEROP_H

class LSOBCCompressibleEulerOp;

#include "lsobcoperator.h"

class LSOBCCompressibleEulerOp : public LSOBCOperator
{
public:

  void access(const LSLayerDataExtrapol &lslde,
              const size_t& field,
              const size_t& abuse_field);

  void operateInner(const LSLayerDataExtrapol& lslde,
                    const size_t& field,
                    const size_t& abuse_field,
                    const real& relax);

  void operateOuter(const LSLayerDataExtrapol& lslde,
                    const size_t& field,
                    const size_t& abuse_field,
                    const real& relax);

};


void LSOBCCompressibleEulerOp::access(const LSLayerDataExtrapol &lslde,
                                      const size_t& i_field,
                                      const size_t& abuse_field)
{
  size_t l_c = lslde.m_Cell;
  //.. pointers
  real* rho   = lslde.m_Data + i_field * lslde.m_FieldSize + 0 * lslde.m_VariableSize;  // 0: rho
  real* rhou  = lslde.m_Data + i_field * lslde.m_FieldSize + 1 * lslde.m_VariableSize;  // 1: rhou
  real* rhov  = lslde.m_Data + i_field * lslde.m_FieldSize + 2 * lslde.m_VariableSize;  // 1: rhov
  real* rhow  = lslde.m_Data + i_field * lslde.m_FieldSize + 3 * lslde.m_VariableSize;  // 1: rhow
  real* rhoE  = lslde.m_Data + i_field * lslde.m_FieldSize + 4 * lslde.m_VariableSize;  // 1: rhoE
  real* rho_acc   = lslde.m_Data + abuse_field * lslde.m_FieldSize + 0 * lslde.m_VariableSize;
  real* rhou_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 1 * lslde.m_VariableSize;
  real* rhov_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 2 * lslde.m_VariableSize;
  real* rhow_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 3 * lslde.m_VariableSize;
  real* rhoE_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 4 * lslde.m_VariableSize;
  //.. access mirror data
  rho_acc[l_c] = 0.;
  rhou_acc[l_c] = 0.;
  rhov_acc[l_c] = 0.;
  rhow_acc[l_c] = 0.;
  rhoE_acc[l_c] = 0.;
  for(size_t i = 0; i < 8; i++) {
    size_t donor = lslde.m_MirrorDonor[i];
    real weight  = lslde.m_MirrorWeight[i];
    rho_acc[l_c]  += weight * rho[donor];
    rhou_acc[l_c] += weight * rhou[donor];
    rhov_acc[l_c] += weight * rhov[donor];
    rhow_acc[l_c] += weight * rhow[donor];
    rhoE_acc[l_c] += weight * rhoE[donor];
  }
}


void LSOBCCompressibleEulerOp::operateInner(const LSLayerDataExtrapol& lslde,
                                            const size_t& i_field,
                                            const size_t& abuse_field,
                                            const real& relax)
{
  size_t l_c = lslde.m_Cell;
  real gxn = lslde.m_Gx;
  real gyn = lslde.m_Gy;
  real gzn = lslde.m_Gz;

  //.. pointers
  real* rho   = lslde.m_Data + i_field * lslde.m_FieldSize + 0 * lslde.m_VariableSize;  // 0: rho
  real* rhou  = lslde.m_Data + i_field * lslde.m_FieldSize + 1 * lslde.m_VariableSize;  // 1: rhou
  real* rhov  = lslde.m_Data + i_field * lslde.m_FieldSize + 2 * lslde.m_VariableSize;  // 1: rhov
  real* rhow  = lslde.m_Data + i_field * lslde.m_FieldSize + 3 * lslde.m_VariableSize;  // 1: rhow
  real* rhoE  = lslde.m_Data + i_field * lslde.m_FieldSize + 4 * lslde.m_VariableSize;  // 1: rhoE
  real* rho_acc   = lslde.m_Data + abuse_field * lslde.m_FieldSize + 0 * lslde.m_VariableSize;
  real* rhou_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 1 * lslde.m_VariableSize;
  real* rhov_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 2 * lslde.m_VariableSize;
  real* rhow_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 3 * lslde.m_VariableSize;
  real* rhoE_acc  = lslde.m_Data + abuse_field * lslde.m_FieldSize + 4 * lslde.m_VariableSize;

  // Apply mirror condition;
  //  * copy rho_acc and rhoE_acc => keep p, T, ... (as also |rho_uvw| is kept)
  //  * keep tangential speed
  //  * reverse normal speed found on access point
  // rhouvw_mirror = rhouvw_acc - 2*(nabla_G * rhouvw) * nabla_G
  //
  //.. Normal projection of rhouvw_acc
  real rho_uvw_n = rhou_acc[l_c] * gxn + rhov_acc[l_c] * gyn + rhow_acc[l_c] * gzn;
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
}


void LSOBCCompressibleEulerOp::operateOuter(const LSLayerDataExtrapol& lslde,
                                            const size_t& i_field,
                                            const size_t& abuse_field,
                                            const real& relax)
{
  // do nothing
}


#endif // LSOBCCOMPRESSIBLEEULEROP_H
