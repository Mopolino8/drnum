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

#include "cartesianpatch.h"
#include "geometrytools.h"

#ifdef WITH_VTK
#include <vtkCellType.h>
#endif


CartesianPatch::CartesianPatch(PatchGrid* patch_grid, size_t num_protectlayers, size_t num_overlaplayers)
  : Patch(patch_grid, num_protectlayers, num_overlaplayers)
{
  m_mypatchtype = 1001;

  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
  //do allways with computeDeltas() m_Interpol_Initialized = false;
  /// @todo need a better eps-handling.
  m_Eps = 1.e-5;
  // NEW NEW_SEEK_EXCEPTION  setNumProtectLayers(m_NumProtectLayers);
}


void CartesianPatch::cellNeighbours (const size_t& l_cell,
                                     vector<size_t>& l_cell_neighbours)
{
  size_t i_cell, j_cell, k_cell;
  ijk(l_cell,
      i_cell, j_cell, k_cell);
  cellNeighbours(i_cell, j_cell, k_cell,
                 l_cell_neighbours);
}


void CartesianPatch::cellNeighbours (const size_t& i_cell,
                                     const size_t& j_cell,
                                     const size_t& k_cell,
                                     vector<size_t>& l_cell_neighbours)
{
  /// @todo if save_index is too slow, can write a faster method, checking bounds

  l_cell_neighbours.clear();
  bool error;
  size_t l_cn;

  int i = int(i_cell);
  int j = int(j_cell);
  int k = int(k_cell);

  l_cn = save_index( i  , j  , k-1,
                    error);
  if (!error) { l_cell_neighbours.push_back(l_cn); };

  l_cn = save_index( i  , j  , k+1,
                    error);
  if (!error) { l_cell_neighbours.push_back(l_cn); };

  l_cn = save_index( i  , j-1, k  ,
                    error);
  if (!error) { l_cell_neighbours.push_back(l_cn); };

  l_cn = save_index( i  , j+1, k  ,
                    error);
  if (!error) { l_cell_neighbours.push_back(l_cn); };

  l_cn = save_index( i-1, j  , k  ,
                    error);
  if (!error) { l_cell_neighbours.push_back(l_cn); };

  l_cn = save_index( i+1, j  , k  ,
                    error);
  if (!error) { l_cell_neighbours.push_back(l_cn); };

}


bool CartesianPatch::readFromFile(istringstream& iss_input, real scale)
{
  bool no_error = Patch::readFromFile(iss_input, scale);
  // number of nodes in block axis directions
  size_t num_i, num_j, num_k;
  iss_input >> num_i;
  iss_input >> num_j;
  iss_input >> num_k;
  // protection exceptions
  size_t numProtXmin, numProtXmax, numProtYmin, numProtYmax, numProtZmin, numProtZmax;
  iss_input >> numProtXmin;
  iss_input >> numProtXmax;
  iss_input >> numProtYmin;
  iss_input >> numProtYmax;
  iss_input >> numProtZmin;
  iss_input >> numProtZmax;
  // physical size of cartesian block
  real ilength, jlength, klength;
  iss_input >> ilength;
  iss_input >> jlength;
  iss_input >> klength;
  // scale length according to IO-scaling factor
  ilength *= m_IOScale;
  jlength *= m_IOScale;
  klength *= m_IOScale;
  // apply patch modifiers
  setSeekExceptions(numProtXmin, numProtXmax, numProtYmin, numProtYmax, numProtZmin, numProtZmax);
  resize(num_i, num_j, num_k);
  buildRegions();  // also in resize
  setupMetrics(ilength, jlength, klength);
  // continue reading solver codes from file
  Patch::readSolverCodes(iss_input);

  /// @todo check before returning "true", see also patch.cpp
  return no_error;
}

//bool CartesianPatch::writeToFile(ifstream &s_mesh)
//{
//  // implementation missing
//}

void CartesianPatch::scaleRefParental(real scfactor)
{
  Patch::scaleRefParental(scfactor);
  m_Lx = scfactor * m_Lx;
  m_Ly = scfactor * m_Ly;
  m_Lz = scfactor * m_Lz;
  computeDeltas();
}


void CartesianPatch::setSeekExceptions(const size_t& num_seekImin, const size_t& num_seekImax,
                                       const size_t& num_seekJmin, const size_t& num_seekJmax,
                                       const size_t& num_seekKmin, const size_t& num_seekKmax)
{
  m_SeekExceptions = true;
  m_NumSeekImin = num_seekImin;
  m_NumSeekImax = num_seekImax;
  m_NumSeekJmin = num_seekJmin;
  m_NumSeekJmax = num_seekJmax;
  m_NumSeekKmin = num_seekKmin;
  m_NumSeekKmax = num_seekKmax;
}


void CartesianPatch::computeDeltas()
{

  m_Dx = m_Lx/m_NumI;
  m_Dy = m_Ly/m_NumJ;
  m_Dz = m_Lz/m_NumK;
  countFlops(3);

  m_InvDx = 1.0/m_Dx;
  m_InvDy = 1.0/m_Dy;
  m_InvDz = 1.0/m_Dz;
  countFlops(3);

  m_xCCMin = 0.5 * m_Dx;
  m_yCCMin = 0.5 * m_Dy;
  m_zCCMin = 0.5 * m_Dz;
  m_xCCMax = m_Lx - 0.5 * m_Dx;
  m_yCCMax = m_Ly - 0.5 * m_Dy;
  m_zCCMax = m_Lz - 0.5 * m_Dz;
  countFlops(6);

  // Set interpolation responsibility box for CC data.
  // Attention: Allow interpolation in total region of cells
  //            in responibility domain (+ 0.) instead of (+0.5).
  //            Reason: Adaptivity. Precision actually handled via epsing.
  // Correction: 2013_06_28: use 0.5 element layer protection for CC-interpolates
  //             to avoid any invalid access causing recursion.
  m_xCCInterMin = m_xCCMin + m_NumProtImin * m_Dx;
  m_yCCInterMin = m_yCCMin + m_NumProtJmin * m_Dy;
  m_zCCInterMin = m_zCCMin + m_NumProtKmin * m_Dz;
  m_xCCInterMax = m_xCCMax - m_NumProtImax * m_Dx;
  m_yCCInterMax = m_yCCMax - m_NumProtJmax * m_Dy;
  m_zCCInterMax = m_zCCMax - m_NumProtKmax * m_Dz;

  // DEBUG: erase below, if above is correct
  m_xCCInterMin = (m_NumProtImin + 0.5) * m_Dx;
  m_yCCInterMin = (m_NumProtJmin + 0.5) * m_Dy;
  m_zCCInterMin = (m_NumProtKmin + 0.5) * m_Dz;
  m_xCCInterMax = m_Lx - (m_NumProtImax + 0.5) * m_Dx;
  m_yCCInterMax = m_Ly - (m_NumProtJmax + 0.5) * m_Dy;
  m_zCCInterMax = m_Lz - (m_NumProtKmax + 0.5) * m_Dz;

  countFlops(6);

  m_EpsDX = m_Dx * m_Eps;
  m_EpsDY = m_Dy * m_Eps;
  m_EpsDZ = m_Dz * m_Eps;
  countFlops(3);
}


void CartesianPatch::setupMetrics(real ilength, real jlength, real klength)
{
  // copy physical size (lengths) of block
  m_Lx = ilength;
  m_Ly = jlength;
  m_Lz = klength;
  // build up increments, etc ...
  computeDeltas();
}


void CartesianPatch::resize(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI  = num_i;
  m_NumJ  = num_j;
  m_NumK  = num_k;
  m_NumJK = num_j*num_k;
  deleteData();
  Patch::resize(m_NumI*m_NumJ*m_NumK);
  buildRegions();
  computeDeltas();
}

void CartesianPatch::buildBoundingBox()
{
  // Test all 8 corners to find max coods in o-system
  //.. prime
  vec3_t v_corner(0, 0, 0);
  vec3_t vo_corner = m_TransformInertial2This.transformReverse(v_corner);
  m_BBoxXYZoMin = vo_corner;
  m_BBoxXYZoMax = vo_corner;

  //.. loop corners
  real xc[2];
  real yc[2];
  real zc[2];
  xc[0] = 0.;
  xc[1] = m_Lx;
  yc[0] = 0.;
  yc[1] = m_Ly;
  zc[0] = 0.;
  zc[1] = m_Lz;

  for (size_t ii=0; ii<2; ii++) {
    v_corner[0] = xc[ii];
    for (size_t jj=0; jj<2; jj++) {
      v_corner[1] = yc[jj];
      for (size_t kk=0; kk<2; kk++) {
        v_corner[2] = zc[kk];
        vo_corner = m_TransformInertial2This.transformReverse(v_corner);
        m_BBoxXYZoMin.minimisePerCoord(vo_corner);
        m_BBoxXYZoMax.maximisePerCoord(vo_corner);
      }
    }
  }
  m_BBoxOk = true;
}


void CartesianPatch::xxyyzzSubCellRaster(const size_t& l_cell, const size_t& lin_mult_res,
                                         vector<vec3_t>& xxyyzz_subcells,
                                         vec3_t& ref_dxxyyzz)
{
  xxyyzz_subcells.clear();

  // i, j, k
  size_t i, j, k;
  ijk(l_cell,
      i, j, k);

  // increment in subcell resolution
  real resolve_factor = 1./float(lin_mult_res);
  real dxx_sub = m_Dx * resolve_factor;
  real dyy_sub = m_Dy * resolve_factor;
  real dzz_sub = m_Dz * resolve_factor;

  // lower corner point in all coords
  real xx_low = i*m_Dx + 0.5 * dxx_sub;
  real yy_low = j*m_Dy + 0.5 * dyy_sub;
  real zz_low = k*m_Dz + 0.5 * dzz_sub;

  // fill subcell coord triples in list
  vec3_t xxyyzz_h;
  for (size_t i_sub = 0; i_sub < lin_mult_res; i_sub ++) {
    xxyyzz_h[0] = xx_low + i_sub * dxx_sub;
    for (size_t j_sub = 0; j_sub < lin_mult_res; j_sub ++) {
      xxyyzz_h[1] = yy_low + j_sub * dyy_sub;
      for (size_t k_sub = 0; k_sub < lin_mult_res; k_sub ++) {
        xxyyzz_h[2] = zz_low + k_sub * dzz_sub;
        xxyyzz_subcells.push_back(xxyyzz_h);
      }
    }
  }
  ref_dxxyyzz[0] = m_Dx;
  ref_dxxyyzz[1] = m_Dy;
  ref_dxxyyzz[2] = m_Dz;
}


void CartesianPatch::buildRegions()
{
  // Example: a 16x8 patch with 2 seek layers and 1 protection layer
  //          with seek exception m_NumSeekImax = m_NumSeekJmin = 0 would
  //          have the following regions below:
  //            S:     seeking zone
  //            S & P: protected (no other patches are allowed to seek here)
  //            D:     donor zone
  //            P & D: core
  //
  //   S S S S S S S S S S S S S S S S
  //   S S S S S S S S S S S S S S S S -- m_JCoreAfterlast
  //   S S P P P P P P P P P P P P P P -- m_DonorZoneAfterlast
  //   S S P D D D D D D D D D D D D D
  //   S S P D D D D D D D D D D D D D
  //   S S P D D D D D D D D D D D D D
  //   S S P D D D D D D D D D D D D D
  //   S S P D D D D D D D D D D D D D
  //   S S P D D D D D D D D D D D D D -- m_JCoreFirst = m_JDonorZoneFirst
  //       | |                         |
  //       | m_IDonorZoneFirst         m_IDonorZoneAfterlast
  //       m_ICoreFirst                m_ICoreAfterlast
  //

  // Seek zone: zone in overlap, seeking data from foreign patches
  if (!m_SeekExceptions) { // no seek exception on whole patch
    m_NumSeekImin = m_NumSeekLayers;
    m_NumSeekImax = m_NumSeekLayers;
    m_NumSeekJmin = m_NumSeekLayers;
    m_NumSeekJmax = m_NumSeekLayers;
    m_NumSeekKmin = m_NumSeekLayers;
    m_NumSeekKmax = m_NumSeekLayers;
  } else {
    // nothing to do: individual seek boundaries have been set on patch reading
  }

  // Protection zone: boundary zone of patch, in which no data access is allowed.
  // Includes seek zone and adds a number of protection layers to the seek zone, if
  // m_NumSeek.... != 0 .
  // If m_NumSeek.... == 0 a boundary is assumed, on which data access is permitted
  // without any protection.
  m_NumProtImin = m_NumSeekImin + m_NumAddProtectLayers;
  m_NumProtImax = m_NumSeekImax + m_NumAddProtectLayers;
  m_NumProtJmin = m_NumSeekJmin + m_NumAddProtectLayers;
  m_NumProtJmax = m_NumSeekJmax + m_NumAddProtectLayers;
  m_NumProtKmin = m_NumSeekKmin + m_NumAddProtectLayers;
  m_NumProtKmax = m_NumSeekKmax + m_NumAddProtectLayers;
  if (m_NumSeekImin == 0) {m_NumProtImin = 0;}
  if (m_NumSeekImax == 0) {m_NumProtImax = 0;}
  if (m_NumSeekJmin == 0) {m_NumProtJmin = 0;}
  if (m_NumSeekJmax == 0) {m_NumProtJmax = 0;}
  if (m_NumSeekKmin == 0) {m_NumProtKmin = 0;}
  if (m_NumSeekKmax == 0) {m_NumProtKmax = 0;}

  // Define Donor region of patch
  // [m_IDonorZoneFirst, m_IDonorZoneAfterlast - 1]
  // [m_JDonorZoneFirst, m_JDonorZoneAfterlast - 1]
  // [m_KDonorZoneFirst, m_KDonorZoneAfterlast - 1]
  m_IDonorZoneFirst = m_NumProtImin;
  m_IDonorZoneAfterlast = m_NumI - m_NumProtImax;
  m_JDonorZoneFirst = m_NumProtJmin;
  m_JDonorZoneAfterlast = m_NumJ - m_NumProtJmax;
  m_KDonorZoneFirst = m_NumProtKmin;
  m_KDonorZoneAfterlast = m_NumK - m_NumProtKmax;

  // Define Core region of patch
  // [m_ICoreFirst, m_ICoreAfterlast - 1]
  // [m_JCoreFirst, m_JCoreAfterlast - 1]
  // [m_KCoreFirst, m_KCoreAfterlast - 1]
  m_ICoreFirst = m_NumSeekImin;
  m_ICoreAfterlast = m_NumI - m_NumSeekImax;
  m_JCoreFirst = m_NumSeekJmin;
  m_JCoreAfterlast = m_NumJ - m_NumSeekJmax;
  m_KCoreFirst = m_NumSeekKmin;
  m_KCoreAfterlast = m_NumK - m_NumSeekKmax;

}


void CartesianPatch::extractSeekCells()
{
  bool any_error = false;
  bool error;
  size_t cell_h;

  /** @todo check!! very prone to cut&paste errors.
    * Can do better than below looping: inserts many duplicates, that will be eliminated later.
    */

  // I-min side
  for (size_t i_cell = 0; i_cell < m_NumSeekImin; i_cell++) {
    for (size_t j_cell = 0; j_cell < m_NumJ; j_cell++) {
      for (size_t k_cell = 0; k_cell < m_NumK; k_cell++) {
        cell_h = save_index(i_cell, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        if (!error) {
          m_ReceiveCells.push_back(cell_h);
        }
      }
    }
  }

  // I-max side
  for (size_t i_cell = m_NumI - m_NumSeekImax; i_cell < m_NumI; i_cell++) {
    for (size_t j_cell = 0; j_cell < m_NumJ; j_cell++) {
      for (size_t k_cell = 0; k_cell < m_NumK; k_cell++) {
        cell_h = save_index(i_cell, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        if (!error) {
          m_ReceiveCells.push_back(cell_h);
        }
      }
    }
  }

  // J-min side
  for (size_t j_cell = 0; j_cell < m_NumSeekJmin; j_cell++) {
    for (size_t k_cell = 0; k_cell < m_NumK; k_cell++) {
      for (size_t i_cell = 0; i_cell < m_NumI; i_cell++) {
        cell_h = save_index(i_cell, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        if (!error) {
          m_ReceiveCells.push_back(cell_h);
        }
      }
    }
  }

  // J-max side
  for (size_t j_cell = m_NumJ - m_NumSeekJmax; j_cell < m_NumJ; j_cell++) {
    for (size_t k_cell = 0; k_cell < m_NumK; k_cell++) {
      for (size_t i_cell = 0; i_cell < m_NumI; i_cell++) {
        cell_h = save_index(i_cell, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        if (!error) {
          m_ReceiveCells.push_back(cell_h);
        }
      }
    }
  }

  // K-min side:
  for (size_t k_cell = 0; k_cell < m_NumSeekKmin; k_cell++) {
    for (size_t i_cell = 0; i_cell < m_NumI; i_cell++) {
      for (size_t j_cell = 0; j_cell < m_NumJ; j_cell++) {
        cell_h = save_index(i_cell, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        if (!error) {
          m_ReceiveCells.push_back(cell_h);
        }
      }
    }
  }

  // K-max side:
  for (size_t k_cell = m_NumK - m_NumSeekKmax; k_cell < m_NumK; k_cell++) {
    for (size_t i_cell = 0; i_cell < m_NumI; i_cell++) {
      for (size_t j_cell = 0; j_cell < m_NumJ; j_cell++) {
        cell_h = save_index(i_cell, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        if (!error) {
          m_ReceiveCells.push_back(cell_h);
        }
      }
    }
  }

  // Eliminate duplicates from m_ReceiveCells
  compactReceiveCellLists();
  m_receiveCells_OK = true;

}


bool CartesianPatch::computeDependencies(const size_t& i_neighbour)
{
  bool found_dependency = false;
  Patch* neighbour_patch = m_neighbours[i_neighbour].first;
  CoordTransformVV trans = m_neighbours[i_neighbour].second;
  // patch bounding normal vectors in coord syst of donor
  vector<vec3_t> nxxyyzz_ijk;
  nxxyyzz_ijk.resize(3);
  nxxyyzz_ijk[0] = trans.transfree(vec3_t(1., 0., 0.));
  nxxyyzz_ijk[1] = trans.transfree(vec3_t(0., 1., 0.));
  nxxyyzz_ijk[2] = trans.transfree(vec3_t(0., 0., 1.));
  // Process all receiving cells and get interpolation sets
  //  - get own coeffs
  //  - transform into system of neighbour patch
  //  - interpolate there if receiving cell hits neighbours core region
  for(size_t ll_rc=0; ll_rc < m_ReceiveCells.size(); ll_rc++) {
    size_t l_rc = m_ReceiveCells[ll_rc];
    vec3_t xyz_rc = xyzCell(l_rc);
    vec3_t xxyyzz_rc = trans.transform(xyz_rc);
    WeightedSet<real> w_set;
    if(m_InterpolateData) {
      //      if(neighbour_patch->computeCCDataInterpolCoeffs(xxyyzz_rc,
      //                                                      w_set)) {
      if(neighbour_patch->computeCCDataInterpolCoeffs_V1(xxyyzz_rc[0], xxyyzz_rc[1], xxyyzz_rc[2],
                                                         w_set)) {
        m_InterCoeffData_WS[i_neighbour].push(ll_rc, l_rc, w_set);  // note: since l_rc are unique, no "add" is required
        m_receive_cell_data_hits[ll_rc]++;
        found_dependency = true;
      }
    }
  }
  return found_dependency;
}


bool CartesianPatch::checkBoxOverlap(const vec3_t& box_xyzo_min, const vec3_t& box_xyzo_max,
                                     const bool& only_core)
{
  bool inside = false;

  // save-side: loop through all cells at a boundary and check inside box
  size_t i_min, i_after_max;
  size_t j_min, j_after_max;
  size_t k_min, k_after_max;

  if (only_core) {
    i_min = m_ICoreFirst;
    j_min = m_JCoreFirst;
    k_min = m_KCoreFirst;
    i_after_max = m_ICoreAfterlast;
    j_after_max = m_JCoreAfterlast;
    k_after_max = m_KCoreAfterlast;
  }

  // loop all cells
  // note: if transformation is too expensive, might also transform box into local
  //       xyz-system and find xyz-coord-limits (all 8 nodes of box!!!)
  for (size_t i_cell = i_min; i_cell < i_after_max; i_cell++) {
    for (size_t j_cell = j_min; j_cell < j_after_max; j_cell++) {
      for (size_t k_cell = k_min; k_cell < k_after_max; k_cell++) {
        vec3_t xyz_cell = xyzCell(i_cell, j_cell, k_cell);
        vec3_t xyzo_cell = m_TransformInertial2This.transformReverse(xyz_cell);
        if (box_xyzo_min[0] <= xyzo_cell[0] &&
            box_xyzo_min[1] <= xyzo_cell[1] &&
            box_xyzo_min[2] <= xyzo_cell[2] &&
            box_xyzo_max[0] >= xyzo_cell[0] &&
            box_xyzo_max[1] >= xyzo_cell[1] &&
            box_xyzo_max[2] >= xyzo_cell[2] ) {
          inside = true;
          return inside;
        }
      }
    }
  }

  return inside;
}


void CartesianPatch::setupInterpolators()
{
  // nothing to be done here.
  /** @note it is likely that setupInterpolators will be meaningful for other
      * patch types. If allways envoqued, shift stuff from computeDeltas back here.
      */
}


bool CartesianPatch::computeCCDataInterpolCoeffs(real x, real y, real z,
                                                 WeightedSet<real>& w_set)
{
  size_t ic_ref, jc_ref, kc_ref;
  bool inside;
  real w_octet[8];

  // Check inside and get reference cell
  inside = xyzToRefInterCell(x, y, z,
                             ic_ref, jc_ref, kc_ref);
  // Get interpol coeff sets
  if(inside) {
    //.. Interpolation weights
    computeCCInterpolWeights(x, y, z,
                             ic_ref, jc_ref, kc_ref,
                             w_octet);
    //.. Fill into w_set
    //.... Set upper cell addresses
    //     NOTE: It is ensured, that ic_ref, jc_ref, kc_ref are in bounds, but
    //           meshes could eventually be flat, (1D, 2D, testing). Incremented
    //           indicees must eventually be shifted back.
    size_t ic_ref_1 = ic_ref + 1;
    if((ic_ref_1 + 1) > m_NumI) {ic_ref_1 = ic_ref;} // must be flat in i-coords
    size_t jc_ref_1 = jc_ref + 1;
    if((jc_ref_1 + 1) > m_NumJ) {jc_ref_1 = jc_ref;} // must be flat in j-coords
    size_t kc_ref_1 = kc_ref + 1;
    if((kc_ref_1 + 1) > m_NumK) {kc_ref_1 = kc_ref;} // must be flat in k-coords
    //.... fill in
    w_set.clearWS();
    w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), w_octet[0]);
    w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref_1), w_octet[1]);
    w_set.pushBack(index(ic_ref,   jc_ref_1, kc_ref  ), w_octet[2]);
    w_set.pushBack(index(ic_ref,   jc_ref_1, kc_ref_1), w_octet[3]);
    w_set.pushBack(index(ic_ref_1, jc_ref,   kc_ref  ), w_octet[4]);
    w_set.pushBack(index(ic_ref_1, jc_ref,   kc_ref_1), w_octet[5]);
    w_set.pushBack(index(ic_ref_1, jc_ref_1, kc_ref  ), w_octet[6]);
    w_set.pushBack(index(ic_ref_1, jc_ref_1, kc_ref_1), w_octet[7]);

    //    w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), w_octet[0]);
    //    w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref+1), w_octet[1]);
    //    w_set.pushBack(index(ic_ref,   jc_ref+1, kc_ref  ), w_octet[2]);
    //    w_set.pushBack(index(ic_ref,   jc_ref+1, kc_ref+1), w_octet[3]);
    //    w_set.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ), w_octet[4]);
    //    w_set.pushBack(index(ic_ref+1, jc_ref,   kc_ref+1), w_octet[5]);
    //    w_set.pushBack(index(ic_ref+1, jc_ref+1, kc_ref  ), w_octet[6]);
    //    w_set.pushBack(index(ic_ref+1, jc_ref+1, kc_ref+1), w_octet[7]);

    //.. Eliminate round of contributors with below eps-weight
    /// @todo m_Eps is used twice, see also computeDeltas. Potentially conflicting.
    w_set.EliminateBelowEps(10*m_Eps, false, false);
    //.. Ensure sum of weights to be 1. to correct eps-errors
    w_set.adjustWeightSumShift(1.);
  }
  return inside;
}

bool CartesianPatch::computeCCGrad1NInterpolCoeffs(real x, real y, real z,
                                                   real nx, real ny, real nz,
                                                   WeightedSet<real>& d_dn)
{
  size_t ic_ref, jc_ref, kc_ref;
  bool inside;
  real w_octet[8];


  /**
    * @todo
    * There is not much advantage to transfer a gradient, if the one computed at the
    * patch of origin requires NumProtectLayers=2 . The total number of overlapping
    * mesh layers would then be 3 (2 on the donor, 1 an the receiving patch). With
    * simple data exchange only, 4 overlapping layers would be required.
    *
    * Optionally, transfer a one sided gradient only. This is possible for NumProtectLayers=1 ,
    * while the receiving patch might reconstruct the central gradient with its own data.
    *
    * NOT YET IMPLEMENTED !!
    *
    */

  //  // May not request a value, if m_NumProtectLayers<2
  //  if(m_NumProtectLayers<2) {
  //    cout << "CartesianPatch::computeCCGrad1NInterpolCoeffs, m_NumProtectLayers = " << m_NumProtectLayers << endl;
  //  }
  // Check inside and get reference cell
  inside = xyzToRefInterCell(x, y, z,
                             ic_ref, jc_ref, kc_ref);
  // Get interpol coeff sets
  if(inside) {
    //.. Interpolation weights
    computeCCInterpolWeights(x, y, z,
                             ic_ref, jc_ref, kc_ref,
                             w_octet);
    //.. Normalize n-vector
    real inv_n_abs = 1./sqrt(nx*nx + ny*ny + nz*nz);
    real nnx = nx * inv_n_abs;
    real nny = ny * inv_n_abs;
    real nnz = nz * inv_n_abs;

    //.. Convenience
    real div_x = nnx/(2*m_Dx);
    real div_y = nny/(2*m_Dy);
    real div_z = nnz/(2*m_Dz);

    d_dn.clearWS();
    //   0:       (ic_ref,   jc_ref,   kc_ref  )
    size_t ii = ic_ref;
    size_t jj = jc_ref;
    size_t kk = kc_ref;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[0]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[0]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[0]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[0]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[0]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[0]*div_z);
    //   1:       (ic_ref,   jc_ref,   kc_ref+1)
    ii = ic_ref;
    jj = jc_ref;
    kk = kc_ref+1;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[1]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[1]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[1]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[1]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[1]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[1]*div_z);
    //   2:       (ic_ref,   jc_ref+1, kc_ref  )
    ii = ic_ref;
    jj = jc_ref+1;
    kk = kc_ref;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[2]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[2]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[2]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[2]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[2]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[2]*div_z);
    //   3:       (ic_ref,   jc_ref+1, kc_ref+1)
    ii = ic_ref;
    jj = jc_ref+1;
    kk = kc_ref+1;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[3]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[3]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[3]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[3]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[3]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[3]*div_z);
    //   4:       (ic_ref+1, jc_ref,   kc_ref  )
    ii = ic_ref+1;
    jj = jc_ref;
    kk = kc_ref;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[4]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[4]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[4]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[4]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[4]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[4]*div_z);
    //   5:       (ic_ref+1, jc_ref,   kc_ref+1)
    ii = ic_ref+1;
    jj = jc_ref;
    kk = kc_ref+1;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[5]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[5]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[5]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[5]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[5]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[5]*div_z);
    //   6:       (ic_ref+1, jc_ref+1, kc_ref  )
    ii = ic_ref+1;
    jj = jc_ref+1;
    kk = kc_ref;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[6]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[6]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[6]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[6]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[6]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[6]*div_z);
    //   7:       (ic_ref+1, jc_ref+1, kc_ref+1)
    ii = ic_ref+1;
    jj = jc_ref+1;
    kk = kc_ref+1;
    d_dn.pushBack(index(ii+1, jj  , kk  ),  w_octet[7]*div_x);
    d_dn.pushBack(index(ii-1, jj  , kk  ), -w_octet[7]*div_x);
    d_dn.pushBack(index(ii  , jj+1, kk  ),  w_octet[7]*div_y);
    d_dn.pushBack(index(ii  , jj-1, kk  ), -w_octet[7]*div_y);
    d_dn.pushBack(index(ii  , jj  , kk+1),  w_octet[7]*div_z);
    d_dn.pushBack(index(ii  , jj  , kk-1), -w_octet[7]*div_z);
    //
    //.. compact to unique addresses
    d_dn.Unify();
    //.. Eliminate round off contributors with below eps-weight
    /// @todo m_Eps is used twice, see also computeDeltas. Potentially conflicting.
    d_dn.EliminateBelowEps(10*m_Eps, true, false);
    //.. Ensure sum of weights to be 0. to correct eps-errors
    /// @todo check function carefully. Sum of Weights must be close to zero anyway.
    d_dn.adjustWeightSumShift(0.);
  }
  return inside;
}


void CartesianPatch::computeCCInterpolWeights(real& x, real& y, real& z,
                                              size_t& ic_ref, size_t& jc_ref, size_t& kc_ref,
                                              real* w_octet)
{
  // Perform a tri-linear interpolation for a point with coordinates (x,y,z)
  // Interpolate values using a box of the inverse mesh build up of geometric centers of cells
  // box corner    cell indicees of box corners
  //   0:         ( ic_ref,   jc_ref,   kc_ref   )
  //   1:         ( ic_ref,   jc_ref,   kc_ref+1 )
  //   2:         ( ic_ref,   jc_ref+1, kc_ref   )
  //   3:         ( ic_ref,   jc_ref+1, kc_ref+1 )
  //   4:         ( ic_ref+1, jc_ref,   kc_ref   )
  //   5:         ( ic_ref+1, jc_ref,   kc_ref+1 )
  //   6:         ( ic_ref+1, jc_ref+1, kc_ref   )
  //   7:         ( ic_ref+1, jc_ref+1, kc_ref+1 )

  // coordinates (local system of this) of reference cell center
  real x_ref = 0.5*m_Dx + ic_ref*m_Dx;
  real y_ref = 0.5*m_Dy + jc_ref*m_Dy;
  real z_ref = 0.5*m_Dz + kc_ref*m_Dz;

  computeInterpolWeights(x, y, z,
                         x_ref, y_ref, z_ref,
                         w_octet);
}

void CartesianPatch::computeNCInterpolWeights(real& x, real& y, real& z,
                                              size_t& in_ref, size_t& jn_ref, size_t& kn_ref,
                                              real* w_octet)
{
  // Perform a tri-linear interpolation for a point with coordinates (x,y,z)
  // Interpolate values using natural mesh nodes
  // box corner    node indicees of box corners
  //   0:         ( in_ref,   jn_ref,   kn_ref   )
  //   1:         ( in_ref,   jn_ref,   kn_ref+1 )
  //   2:         ( in_ref,   jn_ref+1, kn_ref   )
  //   3:         ( in_ref,   jn_ref+1, kn_ref+1 )
  //   4:         ( in_ref+1, jn_ref,   kn_ref   )
  //   5:         ( in_ref+1, jn_ref,   kn_ref+1 )
  //   6:         ( in_ref+1, jn_ref+1, kn_ref   )
  //   7:         ( in_ref+1, jn_ref+1, kn_ref+1 )

  // coordinates (local system of this) of reference node
  real x_ref = in_ref*m_Dx;
  real y_ref = jn_ref*m_Dy;
  real z_ref = kn_ref*m_Dz;

  computeInterpolWeights(x, y, z,
                         x_ref, y_ref, z_ref,
                         w_octet);
}

void CartesianPatch::computeInterpolWeights(real& x, real& y, real& z,
                                            real& x_ref, real& y_ref, real& z_ref,
                                            real* w_octet)
{
  // Perform a tri-linear interpolation for a point with coordinates (x,y,z)
  // in a box given by the frame corner coordinates.

  // NOTE: no range checking. Will perform an false extrapolation, if
  //       (x,y,z) is outside box.

  // box corner:       coords:
  //   0:          ( x_ref,      y_ref,      z_ref      )
  //   1:          ( x_ref,      y_ref,      z_ref+m_Dz )
  //   2:          ( x_ref,      y_ref+m_Dy, z_ref      )
  //   3:          ( x_ref,      y_ref+m_Dy, z_ref+m_Dz )
  //   4:          ( x_ref+m_Dx, y_ref,      z_ref      )
  //   5:          ( x_ref+m_Dx, y_ref,      z_ref+m_Dz )
  //   6:          ( x_ref+m_Dx, y_ref+m_Dy, z_ref      )
  //   7:          ( x_ref+m_Dx, y_ref+m_Dy, z_ref+m_Dz )

  // cartesian distances to box sides
  real low_diff_x_ref = x - x_ref;
  real low_diff_y_ref = y - y_ref;
  real low_diff_z_ref = z - z_ref;
  real up_diff_x_ref = m_Dx - low_diff_x_ref;
  real up_diff_y_ref = m_Dy - low_diff_y_ref;
  real up_diff_z_ref = m_Dz - low_diff_z_ref;

  real inv_cell_vol = m_InvDx * m_InvDy * m_InvDz;

  //.. Init weights
  w_octet[0] = inv_cell_vol;
  w_octet[1] = inv_cell_vol;
  w_octet[2] = inv_cell_vol;
  w_octet[3] = inv_cell_vol;
  w_octet[4] = inv_cell_vol;
  w_octet[5] = inv_cell_vol;
  w_octet[6] = inv_cell_vol;
  w_octet[7] = inv_cell_vol;

  //.. lower x-side, box corners: 0,1,2,3
  w_octet[0] *= up_diff_x_ref;
  w_octet[1] *= up_diff_x_ref;
  w_octet[2] *= up_diff_x_ref;
  w_octet[3] *= up_diff_x_ref;

  //.. upper x-side, box corners: 4,5,6,7
  w_octet[4] *= low_diff_x_ref;
  w_octet[5] *= low_diff_x_ref;
  w_octet[6] *= low_diff_x_ref;
  w_octet[7] *= low_diff_x_ref;

  //.. lower y-side, box corners: 0,1,4,5
  w_octet[0] *= up_diff_y_ref;
  w_octet[1] *= up_diff_y_ref;
  w_octet[4] *= up_diff_y_ref;
  w_octet[5] *= up_diff_y_ref;

  //.. upper y-side, box corners: 2,3,6,7
  w_octet[2] *= low_diff_y_ref;
  w_octet[3] *= low_diff_y_ref;
  w_octet[6] *= low_diff_y_ref;
  w_octet[7] *= low_diff_y_ref;

  //.. lower z-side, box corners: 0,2,4,6
  w_octet[0] *= up_diff_z_ref;
  w_octet[2] *= up_diff_z_ref;
  w_octet[4] *= up_diff_z_ref;
  w_octet[6] *= up_diff_z_ref;

  //.. upper z-side, box corners: 1,3,5,7
  w_octet[1] *= low_diff_z_ref;
  w_octet[3] *= low_diff_z_ref;
  w_octet[5] *= low_diff_z_ref;
  w_octet[7] *= low_diff_z_ref;

}


// void CartesianPatch::computeCCGrad1NInterpolCoeffs(real& x, real& y, real& z,
// 						   size_t& ic_ref, size_t& jc_ref, size_t& kc_ref,
// 						   real nx, real ny, real nz,
// 						   WeightedSet<real>* d_dn)
// {
//   // Perform a tri-linear interpolation for a point with coordinates (x,y,z)
//   // Interpolate values using inverse mesh build up of geometric centers of cells
//   // w_octet       cell
//   //   0:       (ic_ref,   jc_ref,   kc_ref  )
//   //   1:       (ic_ref,   jc_ref,   kc_ref+1)
//   //   2:       (ic_ref,   jc_ref+1, kc_ref  )
//   //   3:       (ic_ref,   jc_ref+1, kc_ref+1)
//   //   4:       (ic_ref+1, jc_ref,   kc_ref  )
//   //   5:       (ic_ref+1, jc_ref,   kc_ref+1)
//   //   6:       (ic_ref+1, jc_ref+1, kc_ref  )
//   //   7:       (ic_ref+1, jc_ref+1, kc_ref+1)
//   //
//   // ATTENTION: Assume inside inverse mesh

//   // Perform a bi-linear interpolation of 4 edge gradients in discrete coord-directions

//   /// @todo Not sure, if this interpolation holds the formal order of a computation

// // X-gradient
// //.. cell-cell difference grad for 0->4
// WeightedSet<real> d_dx_04();
// d_dx_04.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDx);
// d_dx_04.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDx);
// d_dx_04 *= up_diff_y_ref;
// d_dx_04 *= up_diff_z_ref;
// //.. cell-cell difference grad for 1->5
// WeightedSet<real> d_dx_15();
// d_dx_15.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDx);
// d_dx_15.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDx);
// d_dx_15 *= up_diff_y_ref;
// d_dx_15 *= low_diff_z_ref;
// //.. cell-cell difference grad for 2->6
// WeightedSet<real> d_dx_26();
// d_dx_26.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDx);
// d_dx_26.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDx);
// d_dx_26 *= low_diff_y_ref;
// d_dx_26 *= up_diff_z_ref;
// //.. cell-cell difference grad for 3->7
// WeightedSet<real> d_dx_37();
// d_dx_37.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDx);
// d_dx_37.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDx);
// d_dx_37 *= low_diff_y_ref;
// d_dx_37 *= low_diff_z_ref;
// //.. combined x-grad
// WeightedSet<real> d_dx();
// d_dx += d_dx_04;
// d_dx += d_dx_15;
// d_dx += d_dx_26;
// d_dx += d_dx_37;

// // Y-gradient
// //.. cell-cell difference grad for 0->2
// WeightedSet<real> d_dx_02();
// d_dx_02.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDy);
// d_dx_02.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDy);
// d_dx_02 *= up_diff_x_ref;
// d_dx_02 *= up_diff_z_ref;

// // X-gradient
// //.. proj average on "inverse mesh face" i=ic_ref
// WeightedSet<real> face_i_low();
// face_i_low.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), up_diff_y_ref  * up_diff_z_ref ); // 0
// face_i_low.pushBack(index(ic_ref,   jc_ref,   kc_ref+1), up_diff_y_ref  * low_diff_z_ref); // 1
// face_i_low.pushBack(index(ic_ref,   jc_ref+1, kc_ref  ), low_diff_y_ref * up_diff_z_ref ); // 2
// face_i_low.pushBack(index(ic_ref,   jc_ref+1, kc_ref+1), low_diff_y_ref * low_diff_z_ref); // 3
// //.. proj average on "inverse mesh face" i=ic_ref+1
// WeightedSet<real> face_i_up();
// face_i_low.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ), up_diff_y_ref  * up_diff_z_ref ); // 4
// face_i_low.pushBack(index(ic_ref+1, jc_ref,   kc_ref+1), up_diff_y_ref  * low_diff_z_ref); // 5
// face_i_low.pushBack(index(ic_ref+1, jc_ref+1, kc_ref  ), low_diff_y_ref * up_diff_z_ref ); // 6
// face_i_low.pushBack(index(ic_ref+1, jc_ref+1, kc_ref+1), low_diff_y_ref * low_diff_z_ref); // 7
// //.. x-grad as difference from above
// face_i_low *= (m_InvDy * m_InvDz);
// face_i_up  *= (m_InvDy * m_InvDz);
// WeightedSet<real> d_dx();
// d_dx.clearWS();
// d_dx.Concatenate(face_i_up ,  m_InvDx);
// d_dx.Concatenate(face_i_low, -m_InvDx);

// // Y-gradient
// //.. proj average on "inverse mesh face" j=jc_ref
// WeightedSet<real> face_j_low();
// face_i_low.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), up_diff_x_ref  * up_diff_z_ref ); // 0
// face_i_low.pushBack(index(ic_ref,   jc_ref,   kc_ref+1), up_diff_x_ref  * low_diff_z_ref); // 1
// face_i_low.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ), low_diff_x_ref * up_diff_z_ref ); // 4
// face_i_low.pushBack(index(ic_ref+1, jc_ref,   kc_ref+1), low_diff_x_ref * low_diff_z_ref); // 5
// //.. proj average on "inverse mesh face" j=jc_ref+1
// WeightedSet<real> face_j_up();
// face_i_low.pushBack(index(ic_ref,   jc_ref+1, kc_ref  ), up_diff_x_ref  * up_diff_z_ref ); // 2
// face_i_low.pushBack(index(ic_ref,   jc_ref+1, kc_ref+1), up_diff_x_ref  * low_diff_z_ref); // 3
// face_i_low.pushBack(index(ic_ref+1, jc_ref+1, kc_ref  ), low_diff_x_ref * up_diff_z_ref ); // 6
// face_i_low.pushBack(index(ic_ref+1, jc_ref+1, kc_ref+1), low_diff_x_ref * low_diff_z_ref); // 7
// //.. x-grad as difference from above
// face_j_low *= (m_InvDx * m_InvDz);
// face_j_up  *= (m_InvDx * m_InvDz);
// WeightedSet<real> d_dx();
// d_dx.clearWS();
// d_dx.Concatenate(face_j_up ,  m_InvDy);
// d_dx.Concatenate(face_j_low, -m_InvDy);

// // X-gradient
// //.. proj average on "inverse mesh face" i=ic_ref
// WeightedSet<real> face_i_low();
// face_i_low.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), up_diff_y_ref  * up_diff_z_ref ); // 0
// face_i_low.pushBack(index(ic_ref,   jc_ref,   kc_ref+1), up_diff_y_ref  * low_diff_z_ref); // 1
// face_i_low.pushBack(index(ic_ref,   jc_ref+1, kc_ref  ), low_diff_y_ref * up_diff_z_ref ); // 2
// face_i_low.pushBack(index(ic_ref,   jc_ref+1, kc_ref+1), low_diff_y_ref * low_diff_z_ref); // 3
// //.. proj average on "inverse mesh face" i=ic_ref+1
// WeightedSet<real> face_i_up();
// face_i_low.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ), up_diff_y_ref  * up_diff_z_ref ); // 4
// face_i_low.pushBack(index(ic_ref+1, jc_ref,   kc_ref+1), up_diff_y_ref  * low_diff_z_ref); // 5
// face_i_low.pushBack(index(ic_ref+1, jc_ref+1, kc_ref  ), low_diff_y_ref * up_diff_z_ref ); // 6
// face_i_low.pushBack(index(ic_ref+1, jc_ref+1, kc_ref+1), low_diff_y_ref * low_diff_z_ref); // 7
// //.. x-grad as difference from above
// face_i_low *= (m_InvDy * m_InvDz);
// face_i_up  *= (m_InvDy * m_InvDz);
// WeightedSet<real> d_dx();
// d_dx.clearWS();
// d_dx.Concatenate(face_i_up ,  m_InvDx);
// d_dx.Concatenate(face_i_low, -m_InvDx);

// // Directed gradient
// WeightedSet<real> d_
// //.. normalize n-vector
// real inv_n_abs = 1./sqrt(nx*nx + ny*ny + nz*nz);
// nx *= inv_n_abs;
// ny *= inv_n_abs;
// nz *= inv_n_abs;
// //..
// d_dn-clearWS();
// d_dn->Concatenate(d_dx, nx);
// d_dn->Concatenate(d_dy, ny);
// d_dn->Concatenate(d_dz, nz);
// d_dn->Unify();


// new since Sun Jun 30 04:15:41 CEST 2013
bool CartesianPatch::linearCCInterpolCoeffs (const real& s_test,
                                             const size_t& ijk_donor_first, const size_t& ijk_donor_afterlast,
                                             const real& ds, const size_t& num,
                                             size_t& index_cell_ref, real& rs_coeff, int& sector)
{

  // Naming convention:
  //  s:  general coordinate (as given, x,y,z)
  //  rs: relative coordinate (s / ds)

  real rs_test = s_test / ds;

  // Compute limits of s sectors
  real rs0 = float(ijk_donor_first);
  real rs1 = rs0 + 0.5;
  real rs2 = float(ijk_donor_afterlast) - 0.5;
  real rs3 = rs2 + 0.5;

  // Check sector

  if (rs_test < rs0) {
    sector = 0;  // s_test too low
    return false;
  }
  else if (rs_test < rs1) {
    sector = 1;  // s_test in DonorZonePlus (lower side)
  }
  else if (rs_test < rs2) {
    sector = 2;  // s_test in DonorZone
  }
  else if (rs_test < rs3) {
    sector = 3;  // s_test in DonorZonePlus (upper side)
  }
  else {
    sector = 4;  // s_test too large
    return false;
  }

  // Compute reference cell "cell_ref" for interpolation between neighbour
  // cells [cell_ref, cell_ref+1] and linear weight distribution.
  //
  // Decide upon sector
  if (sector == 1) {
    //.. Low side DonorZonePlus
    //   All interpol weight must go to ijk_donor_first, as ijk_donor_first-1 is forbidden.
    //   This corresponds to cells [ijk_donor_first-1, ijk_donor_first] with rs_coeff = 1.
    //   As ijk_donor_first-1 is forbidden, same effect can be obtained shifting cell
    //   to [ijk_donor_first, ijk_donor_first+1] with rs_coeff = 0.
    index_cell_ref = ijk_donor_first;
    rs_coeff = 0.;
  }

  else if (sector == 2) {
    //.. Normal interpolation in DonorZone
    //.... Get reference cell as long int
    long lindex_cell_ref = long(rs_test - 0.5 + float(num)) - num; // num only to ensure rounding down
    //.... Ensure in accessible bounds:
    //     lindex_cell_ref: [ijk_donor_first, (ijk_donor_afterlast - 2)]
    //     gives access:    [ijk_donor_first, (ijk_donor_afterlast - 1)]
    //
    //     Due to precision errors, it is possible that:
    //     * lindex_cell_ref == (ijk_donor_first - 1)
    //       => this would give a weight rs_coeff = ~1, but access to ijk_donor_first - 1 is forbidden
    //     * lindex_cell_ref == (ijk_donor_afterlast - 1)
    //       => would give rs_coeff = ~0, but access (lindex_cell_ref+1) = ijk_donor_afterlast is forbidden
    if (lindex_cell_ref == (ijk_donor_first - 1)) { // might eventually happen due to precision limits
      lindex_cell_ref = ijk_donor_first;
    }
    else if (lindex_cell_ref == ijk_donor_afterlast - 1) { // might eventually happen due to precision limits
      lindex_cell_ref = ijk_donor_afterlast - 2;
    }
    else if (lindex_cell_ref < -1 || lindex_cell_ref > ijk_donor_afterlast - 1) { // must be a mistake
      BUG;
    }
    index_cell_ref = size_t(lindex_cell_ref); // num only to ensure rounding down
    real rs_cell_ref = 0.5 + float(index_cell_ref);  // rs-coordinate of reference cell
    rs_coeff = rs_test - rs_cell_ref;           // weight distribution ( = rs-distance from rs_cell_ref)
  }

  else if (sector == 3) {
    //.. Upper side DonorZonePlus
    //   All weight must go to ijk_donor_afterlast-1, as ijk_donor_afterlast is forbidden.
    //   This corresponds to cells [ijk_donor_afterlast-1, ijk_donor_afterlast] with rs_coeff = 0.
    //   As ijk_donor_afterlast is forbidden, same effect can be obtained shifting cell
    //   to [ijk_donor_afterlast-2, ijk_donor_afterlast-1] with rs_coeff = 1.
    index_cell_ref = ijk_donor_afterlast - 2;
    rs_coeff = 1.;
  }

  else {
    BUG; // never
  }

  // Check result to be sure:
  if (index_cell_ref < ijk_donor_first) {
    BUG;
  }
  if (index_cell_ref > (ijk_donor_afterlast - 2)) {
    BUG;
  }
  if (rs_coeff < (-m_Eps)) {
    BUG;
  }
  if (rs_coeff > (1.+ m_Eps)) {
    BUG;
  }


  return true;
}


// new since Sun Jun 30 04:15:41 CEST 2013
bool CartesianPatch::computeCCDataInterpolCoeffs_V1(real x, real y, real z,
                                                    WeightedSet<real>& w_set)
{
  size_t ic_ref, jc_ref, kc_ref;
  real x_coeff, y_coeff, z_coeff;
  int sector_i, sector_j, sector_k;
  real w_octet[8];

  // Three times linear analysis
  bool inside;

  // I-direction
  inside = linearCCInterpolCoeffs (x,
                                   m_IDonorZoneFirst, m_IDonorZoneAfterlast,
                                   m_Dx, m_NumI,
                                   ic_ref, x_coeff, sector_i);
  if (!inside) {return false;}

  // J-direction
  inside = linearCCInterpolCoeffs (y,
                                   m_JDonorZoneFirst, m_JDonorZoneAfterlast,
                                   m_Dy, m_NumJ,
                                   jc_ref, y_coeff, sector_j);
  if (!inside) {return false;}

  // K-direction
  inside = linearCCInterpolCoeffs (z,
                                   m_KDonorZoneFirst, m_KDonorZoneAfterlast,
                                   m_Dz, m_NumK,
                                   kc_ref, z_coeff, sector_k);
  if (!inside) {return false;}

  // Getting here, means interpolation was OK in all thre dimensions
  // Compute weights for all eight cells potentially contributing
  //
  // Perform a tri-linear interpolation for a point with coordinates (x,y,z)
  // Interpolate values using a box of the inverse mesh build up of geometric centers of cells
  // box corner    cell indicees of box corners
  //   0:         ( ic_ref,   jc_ref,   kc_ref   )
  //   1:         ( ic_ref,   jc_ref,   kc_ref+1 )
  //   2:         ( ic_ref,   jc_ref+1, kc_ref   )
  //   3:         ( ic_ref,   jc_ref+1, kc_ref+1 )
  //   4:         ( ic_ref+1, jc_ref,   kc_ref   )
  //   5:         ( ic_ref+1, jc_ref,   kc_ref+1 )
  //   6:         ( ic_ref+1, jc_ref+1, kc_ref   )
  //   7:         ( ic_ref+1, jc_ref+1, kc_ref+1 )


  w_octet[0] = (1. - x_coeff)  *  (1. - y_coeff)  *  (1. - z_coeff);
  w_octet[1] = (1. - x_coeff)  *  (1. - y_coeff)  *     z_coeff;
  w_octet[2] = (1. - x_coeff)  *     y_coeff      *  (1. - z_coeff);
  w_octet[3] = (1. - x_coeff)  *     y_coeff      *     z_coeff;
  w_octet[4] =    x_coeff      *  (1. - y_coeff)  *  (1. - z_coeff);
  w_octet[5] =    x_coeff      *  (1. - y_coeff)  *     z_coeff;
  w_octet[6] =    x_coeff      *     y_coeff      *  (1. - z_coeff);
  w_octet[7] =    x_coeff      *     y_coeff      *     z_coeff;


  //  //.. Init weights
  //  w_octet[0] = 1.;
  //  w_octet[1] = 1.;
  //  w_octet[2] = 1.;
  //  w_octet[3] = 1.;
  //  w_octet[4] = 1.;
  //  w_octet[5] = 1.;
  //  w_octet[6] = 1.;
  //  w_octet[7] = 1.;

  //  //.. lower x-side, box corners: 0,1,2,3
  //  w_octet[0] *= (1. - x_coeff);
  //  w_octet[1] *= (1. - x_coeff);
  //  w_octet[2] *= (1. - x_coeff);
  //  w_octet[3] *= (1. - x_coeff);

  //  //.. upper x-side, box corners: 4,5,6,7
  //  w_octet[4] *= x_coeff;
  //  w_octet[5] *= x_coeff;
  //  w_octet[6] *= x_coeff;
  //  w_octet[7] *= x_coeff;

  //  //.. lower y-side, box corners: 0,1,4,5
  //  w_octet[0] *= (1. - y_coeff);
  //  w_octet[1] *= (1. - y_coeff);
  //  w_octet[4] *= (1. - y_coeff);
  //  w_octet[5] *= (1. - y_coeff);

  //  //.. upper y-side, box corners: 2,3,6,7
  //  w_octet[2] *= y_coeff;
  //  w_octet[3] *= y_coeff;
  //  w_octet[6] *= y_coeff;
  //  w_octet[7] *= y_coeff;

  //  //.. lower z-side, box corners: 0,2,4,6
  //  w_octet[0] *= (1. - z_coeff);
  //  w_octet[2] *= (1. - z_coeff);
  //  w_octet[4] *= (1. - z_coeff);
  //  w_octet[6] *= (1. - z_coeff);

  //  //.. upper z-side, box corners: 1,3,5,7
  //  w_octet[1] *= z_coeff;
  //  w_octet[3] *= z_coeff;
  //  w_octet[5] *= z_coeff;
  //  w_octet[7] *= z_coeff;

  // Fill into w_set
  //.. Set upper cell addresses
  //   NOTE: It is ensured, that ic_ref, jc_ref, kc_ref are in bounds, but
  //         meshes could eventually be flat, (1D, 2D, testing). Incremented
  //         indicees must eventually be shifted back.
  size_t ic_ref_1 = ic_ref + 1;
  if((ic_ref_1 + 1) > m_NumI) {ic_ref_1 = ic_ref;} // must be flat in i-coords
  size_t jc_ref_1 = jc_ref + 1;
  if((jc_ref_1 + 1) > m_NumJ) {jc_ref_1 = jc_ref;} // must be flat in j-coords
  size_t kc_ref_1 = kc_ref + 1;
  if((kc_ref_1 + 1) > m_NumK) {kc_ref_1 = kc_ref;} // must be flat in k-coords
  //.... fill in
  w_set.clearWS();
  w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), w_octet[0]);
  w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref_1), w_octet[1]);
  w_set.pushBack(index(ic_ref,   jc_ref_1, kc_ref  ), w_octet[2]);
  w_set.pushBack(index(ic_ref,   jc_ref_1, kc_ref_1), w_octet[3]);
  w_set.pushBack(index(ic_ref_1, jc_ref,   kc_ref  ), w_octet[4]);
  w_set.pushBack(index(ic_ref_1, jc_ref,   kc_ref_1), w_octet[5]);
  w_set.pushBack(index(ic_ref_1, jc_ref_1, kc_ref  ), w_octet[6]);
  w_set.pushBack(index(ic_ref_1, jc_ref_1, kc_ref_1), w_octet[7]);

  //.. Eliminate round of contributors with below eps-weight
  /// @todo m_Eps is used twice, see also computeDeltas. Potentially conflicting.
  w_set.EliminateBelowEps(10*m_Eps, false, false);
  //.. Ensure sum of weights to be 1. to correct eps-errors
  w_set.adjustWeightSumShift(1.);

  return true;
}


bool CartesianPatch::xyzToRefInterCell(real x, real y, real z,
                                       size_t& ic_ref, size_t& jc_ref, size_t& kc_ref)
{
  bool inside = true;
  //
  // Check limits. Allow x,y,z in eps-bounds of interpolation range as well.
  // NOTE: if protective layers m_numProt.... == 0 (possible as boundary exception),
  // locations of (x,y,z) outside the box bounding CCs are possible.
  //
  if(x < m_xCCInterMin-m_EpsDX) {ic_ref = 0; inside = false;};
  if(y < m_yCCInterMin-m_EpsDY) {jc_ref = 0; inside = false;};
  if(z < m_zCCInterMin-m_EpsDZ) {kc_ref = 0; inside = false;};
  if(x > m_xCCInterMax+m_EpsDX) {ic_ref = m_NumI-1; inside = false;};
  if(y > m_yCCInterMax+m_EpsDY) {jc_ref = m_NumJ-1; inside = false;};
  if(z > m_zCCInterMax+m_EpsDZ) {kc_ref = m_NumK-1; inside = false;};
  //
  // Indicees pick, if inside
  if(inside) {
    //..Position relative to lowest CC-coords (indicees hashing)
    real ric_ref = (x - m_xCCMin) / m_Dx;
    real rjc_ref = (y - m_yCCMin) / m_Dy;
    real rkc_ref = (z - m_zCCMin) / m_Dz;
    //.. Correct relative position to ensure indicees are in range
    //   (lower bounds correction: ensure non-negative)
    if(ric_ref < m_Eps) {ric_ref = m_Eps;}
    if(rjc_ref < m_Eps) {rjc_ref = m_Eps;}
    if(rkc_ref < m_Eps) {rkc_ref = m_Eps;}
    //.. Convert to size_t
    ic_ref = size_t(ric_ref);
    jc_ref = size_t(rjc_ref);
    kc_ref = size_t(rkc_ref);
    //
    //.. Correct reference cell indices to ensure inside donor region
    //   NOTE: Later interpolation will be between
    //          [ic_ref , (ic_ref + 1)]
    //          [jc_ref , (jc_ref + 1)]
    //          [kc_ref , (kc_ref + 1)]
    //   => ic_ref, jc_ref, kc_ref restricted to
    //          [m_IDonorZoneFirst , m_IDonorZoneAfterlast-2]
    //          [m_JDonorZoneFirst , m_JDonorZoneAfterlast-2]
    //          [m_KDonorZoneFirst , m_KDonorZoneAfterlast-2]
    if(ic_ref < m_IDonorZoneFirst) {
      ic_ref = m_IDonorZoneFirst;}
    if(jc_ref < m_JDonorZoneFirst) {
      jc_ref = m_JDonorZoneFirst;}
    if(kc_ref < m_KDonorZoneFirst) {
      kc_ref = m_KDonorZoneFirst;}
    if(ic_ref > m_IDonorZoneAfterlast-2) {
      ic_ref = m_IDonorZoneAfterlast-2;}
    if(jc_ref > m_JDonorZoneAfterlast-2) {
      jc_ref = m_JDonorZoneAfterlast-2;}
    if(kc_ref > m_KDonorZoneAfterlast-2) {
      kc_ref = m_KDonorZoneAfterlast-2;}
  }


  return inside;
}


bool CartesianPatch::xyzToRefCell(real x, real y, real z,
                                  size_t& ic_ref, size_t& jc_ref, size_t& kc_ref)
{
  bool inside = true;

  /// @todo Eps-stuff is quite ugly down here

  // Check inside limits. Allow x,y,z in eps-bounds as well.
  if(x < m_xCCMin-m_EpsDX) {ic_ref = 0; inside = false;};
  if(y < m_yCCMin-m_EpsDY) {jc_ref = 0; inside = false;};
  if(z < m_zCCMin-m_EpsDZ) {kc_ref = 0; inside = false;};
  if(x > m_xCCMax+m_EpsDX) {ic_ref = m_NumI-1; inside = false;};
  if(y > m_yCCMax+m_EpsDY) {jc_ref = m_NumJ-1; inside = false;};
  if(z > m_zCCMax+m_EpsDZ) {kc_ref = m_NumK-1; inside = false;};
  // get reference cell, if inside
  if(inside) {
    //.. Manipulate x,y,z if these coords depass the borders
    //   NOTE: inside==true means that x,y,z is in eps-bounds anyway
    //         if close, shift x,y,z to inside to ensure correct address pick
    if(x < m_xCCMin+m_EpsDX) {x = m_xCCMin+m_EpsDX;}
    if(y < m_yCCMin+m_EpsDY) {y = m_yCCMin+m_EpsDY;}
    if(z < m_zCCMin+m_EpsDZ) {z = m_zCCMin+m_EpsDZ;}
    if(x > m_xCCMax-m_EpsDX) {x = m_xCCMax-m_EpsDX;}
    if(y > m_yCCMax-m_EpsDY) {y = m_yCCMax-m_EpsDY;}
    if(z > m_zCCMax-m_EpsDZ) {z = m_zCCMax-m_EpsDZ;}
    //.. Indicees pick
    //   NOTE: Getting here, means "inside" anyway.
    //.... Position relative to lowest CC-coords (indicees hashing)
    real ric_ref = (x - m_xCCMin) / m_Dx;
    real rjc_ref = (y - m_yCCMin) / m_Dy;
    real rkc_ref = (z - m_zCCMin) / m_Dz;
    //.... Correct relative position to ensure indicees are in range
    //     (lower bounds correction: ensure non-negative)
    if(ric_ref < m_Eps) {ric_ref = m_Eps;}
    if(rjc_ref < m_Eps) {rjc_ref = m_Eps;}
    if(rkc_ref < m_Eps) {rkc_ref = m_Eps;}
    //.... Convert to size_t
    ic_ref = size_t(ric_ref);
    jc_ref = size_t(rjc_ref);
    kc_ref = size_t(rkc_ref);
    //.... Correct to upper bounds index limits
    if(ic_ref > m_NumI-2) {ic_ref = m_NumI-2;}
    if(jc_ref > m_NumJ-2) {jc_ref = m_NumJ-2;}
    if(kc_ref > m_NumK-2) {kc_ref = m_NumK-2;}
  }
  return inside;
}

bool CartesianPatch::xyzToRefNode(real x, real y, real z,
                                  size_t& in_ref, size_t& jn_ref, size_t& kn_ref)
{
  bool inside = true;

  /// @todo Eps-stuff is quite ugly down here

  // check limits
  if(x < -m_EpsDX) {in_ref = 0; inside = false;};
  if(y < -m_EpsDY) {jn_ref = 0; inside = false;};
  if(z < -m_EpsDZ) {kn_ref = 0; inside = false;};
  if(x > m_Lx+m_EpsDX) {in_ref = m_NumI-1; inside = false;};
  if(y > m_Ly+m_EpsDY) {jn_ref = m_NumJ-1; inside = false;};
  if(z > m_Lz+m_EpsDZ) {kn_ref = m_NumK-1; inside = false;};

  // get reference cell, if inside
  if(inside) {
    //.. Manipulate x,y,z if these coords depass the borders
    //   NOTE: inside==true means that x,y,z is in eps-bounds anyway
    //         if close, shift x,y,z to inside to ensure correct address pick
    if(x < m_EpsDX) {x = m_EpsDX;}
    if(y < m_EpsDY) {y = m_EpsDY;}
    if(z < m_EpsDZ) {z = m_EpsDZ;}
    if(x > m_Lx-m_EpsDX) {x = m_Lx-m_EpsDX;}
    if(y > m_Ly-m_EpsDY) {y = m_Ly-m_EpsDY;}
    if(z > m_Lz-m_EpsDZ) {z = m_Lz-m_EpsDZ;}
    //.. Position relative to lowest CC-coords.
    real rin_ref = (x - m_xCCMin) / m_Dx;
    real rjn_ref = (y - m_yCCMin) / m_Dy;
    real rkn_ref = (z - m_zCCMin) / m_Dz;
    //    iin_ref = lrint(floor(rin_ref));
    //    jjn_ref = lrint(floor(rjn_ref));
    //    kkn_ref = lrint(floor(rkn_ref));
    //.. Compute reference address. Ensure positive.
    ///@todo: Check, if this is aceptable from robustnes
    //   NOTE: Due to Eps-errors, negative values for rin_ref, rjn_ref, rkn_ref are possible
    //         if Eps_DXYZ==0.
    in_ref = size_t(rin_ref);
    jn_ref = size_t(rjn_ref);
    kn_ref = size_t(rkn_ref);
    //.. Error checking, actually impossible: Be sure to be in addressing range
    bool error = false;
    if(rin_ref < 0.) {
      cout << "CartesianPatch::xyzToRefNode, i_min"; error = true;
      in_ref = 0;
    }
    if(rjn_ref < 0.) {
      cout << "CartesianPatch::xyzToRefNode, j_min"; error = true;
      jn_ref = 0;
    }
    if(rkn_ref < 0.) {
      cout << "CartesianPatch::xyzToRefNode, k_min"; error = true;
      kn_ref = 0;
    }
    if(in_ref > m_NumI-2) {
      cout << "CartesianPatch::xyzToRefNode, i_max"; error = true;
    }
    if(jn_ref > m_NumJ-2) {
      cout << "CartesianPatch::xyzToRefNode, j_max"; error = true;
    }
    if(kn_ref > m_NumK-2) {
      cout << "CartesianPatch::xyzToRefNode, k_max"; error = true;
    }
    if(error) {
      /// @todo may put this in debug condition
      cout << endl;
      cout << " m_xyzCCMin = " << m_xCCMin << ", " << m_yCCMin << ", " << m_zCCMin << endl;
      cout << " m_xyzCCMax = " << m_xCCMax << ", " << m_yCCMax << ", " << m_zCCMax << endl;
      cout << " (x,y,z)    = " << x        << ", " << y        << ", " << z             << endl;
      exit(EXIT_FAILURE);
    }
  }
  return inside;
}


real CartesianPatch::computeMinChLength()
{
  real min_ch_len = m_Dx;
  if(min_ch_len > m_Dy) {
    min_ch_len = m_Dy;
  }
  if(min_ch_len > m_Dz) {
    min_ch_len = m_Dz;
  }
  return min_ch_len;
}


void CartesianPatch::cellOverFaceNeighbours (const size_t& l_cell,
                                             vector<size_t>& l_cell_neighbours)
{
  size_t i_cell, j_cell, k_cell;

  ijk(l_cell,
      i_cell, j_cell, k_cell);

  cellOverFaceNeighbours(i_cell, j_cell, k_cell,
                         l_cell_neighbours);
}


void CartesianPatch::cellOverFaceNeighbours (const size_t& i,
                                             const size_t& j,
                                             const size_t& k,
                                             vector<size_t>& l_cell_neighbours)
{

  l_cell_neighbours.clear();

  // check face neighbours
  // NOTE: size_t(0-1) is huge and will be detected by save_index
  bool error;
  size_t l_fn;

  l_fn = save_index(i  , j  , k-1,
                    error);
  if (!error) l_cell_neighbours.push_back(l_fn);

  l_fn = save_index(i  , j  , k+1,
                    error);
  if (!error) l_cell_neighbours.push_back(l_fn);

  l_fn = save_index(i  , j-1, k  ,
                    error);
  if (!error) l_cell_neighbours.push_back(l_fn);

  l_fn = save_index(i  , j+1, k  ,
                    error);
  if (!error) l_cell_neighbours.push_back(l_fn);

  l_fn = save_index(i-1, j  , k  ,
                    error);
  if (!error) l_cell_neighbours.push_back(l_fn);

  l_fn = save_index(i+1, j  , k  ,
                    error);
  if (!error) l_cell_neighbours.push_back(l_fn);

}


void CartesianPatch::cellOverNodeNeighbours (const size_t& l_cell,
                                             vector<size_t>& l_cell_neighbours)
{
  size_t i_c, j_c, k_c;

  ijk(l_cell,
      i_c, j_c, k_c);

  cellOverNodeNeighbours(i_c, j_c, k_c,
                         l_cell_neighbours);
}


void CartesianPatch::cellOverNodeNeighbours (const size_t& i_c,
                                             const size_t& j_c,
                                             const size_t& k_c,
                                             vector<size_t>& l_cell_neighbours)
{
  l_cell_neighbours.clear();

  // check face neighbours
  bool error;
  size_t l_fn;

  for (int i_shift = -1; i_shift <= 1; i_shift++) {
    for (int j_shift = -1; j_shift <= 1; j_shift++) {
      for (int k_shift = -1; k_shift <= 1; k_shift++) {
        if (i_shift!=0 || j_shift!= 0 || k_shift!=0) {  // any non-0, forbids 0,0,0
          l_fn = save_index(i_c + i_shift, j_c + j_shift, k_c + k_shift,
                            error);
          if(!error) {
            l_cell_neighbours.push_back(l_fn);
          }
        }
      }
    }
  }
}


void CartesianPatch::computeNablaVar(const size_t& field_index, const size_t& var_index, const size_t& l_cell,
                                     real& dvar_dx, real& dvar_dy, real& dvar_dz)
{
  size_t i, j, k;
  size_t l_plus, l_minus;
  bool skip_plus, skip_minus;

  real* var = getVariable(field_index, var_index);

  ijk(l_cell,
      i, j, k);

  // dvar_dx
  l_plus  = save_index(i+1, j  , k  ,
                       skip_plus);
  l_minus = save_index(i-1, j  , k  ,
                       skip_minus);
  if (!skip_plus && !skip_minus) {  // go central
    dvar_dx = (var[l_plus] - var[l_minus]) / (2.*m_Dx);
  }
  else if (skip_plus) {  // one sided i-1, i
    dvar_dx = (var[l_cell] - var[l_minus]) / m_Dx;
  }
  else if (skip_minus) {  // one sided i, i+1
    dvar_dx = (var[l_plus] - var[l_cell]) / m_Dx;
  }
  else {
    BUG;
  }

  // dvar_dy
  l_plus  = save_index(i  , j+1, k  ,
                       skip_plus);
  l_minus = save_index(i  , j-1, k  ,
                       skip_minus);
  if (!skip_plus && !skip_minus) {  // go central
    dvar_dy = (var[l_plus] - var[l_minus]) / (2.*m_Dy);
  }
  else if (skip_plus) {  // one sided j-1, j
    dvar_dy = (var[l_cell] - var[l_minus]) / m_Dy;
  }
  else if (skip_minus) {  // one sided j, j+1
    dvar_dy = (var[l_plus] - var[l_cell]) / m_Dy;
  }
  else {
    BUG;
  }

  // dvar_dz
  l_plus  = save_index(i  , j  , k+1,
                       skip_plus);
  l_minus = save_index(i  , j  , k-1,
                       skip_minus);
  if (!skip_plus && !skip_minus) {  // go central
    dvar_dz = (var[l_plus] - var[l_minus]) / (2.*m_Dz);
  }
  else if (skip_plus) {  // one sided k-1, k
    dvar_dz = (var[l_cell] - var[l_minus]) / m_Dz;
  }
  else if (skip_minus) {  // one sided k, k+1
    dvar_dz = (var[l_plus] - var[l_cell]) / m_Dz;
  }
  else {
    BUG;
  }
}


#ifdef WITH_VTK
vtkSmartPointer<vtkDataSet> CartesianPatch::createVtkDataSet(size_t i_field, const PostProcessingVariables &proc_vars)
{
  /// @todo: must transform output for non o-aligned patches !!!
  /// @todo: does vtk know something like a general transformation in space ???
  // Transform: use only linear transformations at present

  vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();

  size_t i_start = m_NumSeekImin;
  size_t i_stop  = m_NumI - m_NumSeekImax;
  size_t j_start = m_NumSeekJmin;
  size_t j_stop  = m_NumJ - m_NumSeekJmax;
  size_t k_start = m_NumSeekKmin;
  size_t k_stop  = m_NumK - m_NumSeekKmax;

  grid->SetDimensions(i_stop - i_start + 1, j_stop - j_start + 1, k_stop - k_start + 1);
  size_t num_tuples = (i_stop - i_start)*(j_stop - j_start)*(k_stop - k_start);

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  for (size_t k = k_start; k <= k_stop; ++k) {
    for (size_t j = j_start; j <= j_stop; ++j) {
      for (size_t i = i_start; i <= i_stop; ++i) {
        vec3_t xyz_p;
        vec3_t xyzo_p;
        xyz_p[0] = i * m_Dx;
        xyz_p[1] = j * m_Dy;
        xyz_p[2] = k * m_Dz;
        xyzo_p = m_TransformInertial2This.transformReverse(xyz_p);
        points->InsertNextPoint(xyzo_p.data());
      }
    }
  }

  grid->SetPoints(points);

  real* raw_var = new real [numVariables()];
  for (int i_var = 0; i_var < proc_vars.numScalars(); ++i_var) {
    vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
    var->SetName(proc_vars.getScalarName(i_var).c_str());
    var->SetNumberOfValues(num_tuples);
    grid->GetCellData()->AddArray(var);
    vtkIdType id = 0;
    for (size_t k = k_start; k < k_stop; ++k) {
      for (size_t j = j_start; j < j_stop; ++j) {
        for (size_t i = i_start; i < i_stop; ++i) {
          getVarDim(numVariables(), i_field, i, j, k, raw_var);
          var->SetValue(id, proc_vars.getScalar(i_var, raw_var));
          ++id;
        }
      }
    }
  }
  for (int i_var = 0; i_var < proc_vars.numVectors(); ++i_var) {
    vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
    var->SetName(proc_vars.getVectorName(i_var).c_str());
    var->SetNumberOfComponents(3);
    var->SetNumberOfTuples(num_tuples);
    grid->GetCellData()->AddArray(var);
    vtkIdType id = 0;
    for (size_t k = k_start; k < k_stop; ++k) {
      for (size_t j = j_start; j < j_stop; ++j) {
        for (size_t i = i_start; i < i_stop; ++i) {
          getVarDim(numVariables(), i_field, i, j, k, raw_var);
          vec3_t v = proc_vars.getVector(i_var, raw_var);
          v = m_TransformInertial2This.transfreeReverse(v);
          float vf[3];
          vf[0] = v[0]; vf[1] = v[1]; vf[2] = v[2];
          var->SetTuple(id, vf);
          ++id;
        }
      }
    }
  }
  delete [] raw_var;
  return grid;
}

vtkSmartPointer<vtkUnstructuredGrid> CartesianPatch::createVtkGridForCells(const list<size_t> &cells)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  size_t num_nodes = 0;
  vector<int> n2n((sizeI() + 1)*(sizeJ() + 1)*(sizeK() + 1), -1);
  vector<vec3_t> x_node((sizeI() + 1)*(sizeJ() + 1)*(sizeK() + 1));
  for (list<size_t>::const_iterator i_cells = cells.begin(); i_cells != cells.end(); ++i_cells) {
    size_t i, j, k;
    ijk(*i_cells, i, j, k);
    vector<size3_t> N(8);
    N[0] = size3_t(i,   j,   k  );
    N[1] = size3_t(i+1, j,   k  );
    N[2] = size3_t(i+1, j+1, k  );
    N[3] = size3_t(i,   j+1, k  );
    N[4] = size3_t(i,   j,   k+1);
    N[5] = size3_t(i+1, j,   k+1);
    N[6] = size3_t(i+1, j+1, k+1);
    N[7] = size3_t(i,   j+1, k+1);
    for (size_t i_node = 0; i_node < 8; ++i_node) {
      size_t id_node = nodeIndex(N[i_node]);
      if (n2n[id_node] == -1) {
        n2n[id_node] = num_nodes;
        x_node[id_node][0] = N[i_node].i * m_Dx;
        x_node[id_node][1] = N[i_node].j * m_Dy;
        x_node[id_node][2] = N[i_node].k * m_Dz;
        x_node[id_node] = m_TransformInertial2This.transformReverse(x_node[id_node]);
        ++num_nodes;
      }
    }
  }
  points->SetNumberOfPoints(num_nodes);
  grid->SetPoints(points);
  grid->Allocate(cells.size());
  for (list<size_t>::const_iterator i_cells = cells.begin(); i_cells != cells.end(); ++i_cells) {
    size_t i, j, k;
    ijk(*i_cells, i, j, k);
    vector<size_t> id_node(8);
    id_node[0] = nodeIndex(i,   j,   k  );
    id_node[1] = nodeIndex(i+1, j,   k  );
    id_node[2] = nodeIndex(i+1, j+1, k  );
    id_node[3] = nodeIndex(i,   j+1, k  );
    id_node[4] = nodeIndex(i,   j,   k+1);
    id_node[5] = nodeIndex(i+1, j,   k+1);
    id_node[6] = nodeIndex(i+1, j+1, k+1);
    id_node[7] = nodeIndex(i,   j+1, k+1);
    vtkIdType num_pts = 8, pts[8];
    for (size_t i_node = 0; i_node < 8; ++i_node) {
      grid->GetPoints()->SetPoint(n2n[id_node[i_node]], x_node[id_node[i_node]].data());
      pts[i_node] = n2n[id_node[i_node]];
    }
    grid->InsertNextCell(VTK_HEXAHEDRON, num_pts, pts);
  }
  return grid;
}

#endif



void CartesianPatch::writeData(QString base_data_filename, int count)
{
  QString str_filename = base_data_filename;
  if (count >= 0) {
    QString str_myindex;
    QString str_num;
    str_myindex.setNum(m_MyIndex);
    str_num.setNum(count);
    while (str_myindex.size() < 6) {
      str_myindex = "0" + str_myindex;
    }
    while (str_num.size() < 6) {
      str_num = "0" + str_num;
    }
    str_filename = base_data_filename += "_ip" + str_myindex + str_num;
  }

  /// @todo Need better data output construction.
  //  cout << "writing to file \"" << qPrintable(str_filename) << ".vtr\"" << endl;
  //  CompressibleVariables<PerfectGas> cvars;
  //  patch.writeToVtk(0, str_filename, cvars);
  BUG;
}

int CartesianPatch::findCell(vec3_t xo)
{
  vec3_t x = m_TransformInertial2This.transform(xo);
  vec3_t x0(0, 0, 0);
  vec3_t x1(sizeI()*dx(), sizeJ()*dy(), sizeK()*dz());
  if (!GeometryTools::isInsideCartesianBox(x, x0, x1)) {
    return -1;
  }
  int i = x[0]/dx();
  int j = x[1]/dy();
  int k = x[2]/dz();
  return index(i, j, k);
}

list<size_t> CartesianPatch::getNeighbours(size_t idx)
{
  list<size_t> neigh;
  size_t i, j, k;
  ijk(idx, i, j, k);
  if (checkRange(i-1, j, k)) {
    neigh.push_back(index(i-1, j, k));
  }
  if (checkRange(i+1, j, k)) {
    neigh.push_back(index(i+1, j, k));
  }
  if (checkRange(i, j-1, k)) {
    neigh.push_back(index(i, j-1, k));
  }
  if (checkRange(i, j+1, k)) {
    neigh.push_back(index(i, j+1, k));
  }
  if (checkRange(i, j, k-1)) {
    neigh.push_back(index(i, j, k-1));
  }
  if (checkRange(i, j, k+1)) {
    neigh.push_back(index(i, j, k+1));
  }
  return neigh;
}
