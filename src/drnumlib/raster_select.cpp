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
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "structuredhexraster.h"

StructuredHexRaster::StructuredHexRaster()
{
  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
  //m_CellLink = new List(1,1);
  //m_CellLink = NULL;
  m_Eps = 1.e-6;  /// @todo need a better eps-handling.
}


void StructuredHexRaster::computeDeltas()
{
  m_Lx = sqrt(sqr(m_Uxo) + sqr(m_Uyo) + sqr(m_Uzo));
  m_Ly = sqrt(sqr(m_Vxo) + sqr(m_Vyo) + sqr(m_Vzo));
  m_Lz = sqrt(sqr(m_Wxo) + sqr(m_Wyo) + sqr(m_Wzo));
  countFlops(15);
  countSqrts(3);

  m_Dx = m_Lx/m_NumI;
  m_Dy = m_Ly/m_NumJ;
  m_Dz = m_Lz/m_NumK;
  countFlops(3);

  m_InvDx = 1.0/m_Dx;
  m_InvDy = 1.0/m_Dx;
  m_InvDz = 1.0/m_Dx;
  countFlops(3);

  m_xCCMin = 0.5 * m_Dx;
  m_yCCMin = 0.5 * m_Dy;
  m_zCCMin = 0.5 * m_Dz;
  m_xCCMax = m_Lx - 0.5 * m_Dx;
  m_yCCMax = m_Ly - 0.5 * m_Dy;
  m_zCCMax = m_Lz - 0.5 * m_Dz;
  countFlops(6);

  m_xCCInterMin = (m_NumProtectLayers + 0.5) * m_Dx;
  m_yCCInterMin = (m_NumProtectLayers + 0.5) * m_Dy;
  m_zCCInterMin = (m_NumProtectLayers + 0.5) * m_Dz;
  m_xCCInterMax = m_Lx - (m_NumProtectLayers + 0.5) * m_Dx;
  m_yCCInterMax = m_Ly - (m_NumProtectLayers + 0.5) * m_Dy;
  m_zCCInterMax = m_Lz - (m_NumProtectLayers + 0.5) * m_Dz;
  countFlops(6);

  m_EpsDX = m_Lx * m_Eps;
  m_EpsDY = m_Ly * m_Eps;
  m_EpsDZ = m_Lz * m_Eps;
  countFlops(3);
}

void StructuredHexRaster::setupAligned(real xo1, real yo1, real zo1, real xo2, real yo2, real zo2)
{

  /** @todo Use a more general alignment method with trafo-matrix and internal coord system.
   *  I assume these are co-aligned here. */

  m_Xo = xo1;
  m_Yo = yo1;
  m_Zo = zo1;

  m_Uxo = xo2 - xo1;
  m_Uyo = 0;
  m_Uzo = 0;
  countFlops(1);

  m_Vxo = 0;
  m_Vyo = yo2 - yo1;
  m_Vzo = 0;
  countFlops(1);

  m_Wxo = 0;
  m_Wyo = 0;
  m_Wzo = zo2 - zo1;
  countFlops(1);

  // position relative to inertial system
  vec3_t shift_inertial2this;
  shift_inertial2this[0] = -m_Xo;
  shift_inertial2this[1] = -m_Yo;
  shift_inertial2this[2] = -m_Zo;
  m_transformInertial2This.setVector(shift_inertial2this);
/// @todo rotation NOT implemented!!!

  computeDeltas();

//  Transformation t;    /// @todo keep for compatibility
//  t.setVector(vec3_t(xo1, yo1, zo1));
//  setTransformation(t.inverse());
}


void StructuredHexRaster::resize(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  computeDeltas();
//  if(m_CellLink) {
//    delete m_CellLink;
//  }
//  size_t num_cells = m_NumI * m_NumJ * m_NumK;
//  m_CellLink = new List(num_cells, num_cells/10);  /// @todo Our version of "List" sets a minimum increment
}


void StructuredHexRaster::setNumProtectLayers(size_t num_protectlayers)
{
  m_NumProtectLayers = num_protectlayers;
}


void StructuredHexRaster::buildBoundingBox()
{
  vec3_t bbox_xyz_min(0., 0., 0.);
  vec3_t bbox_xyz_max(m_Lx, m_Ly, m_Lz);
  m_bbox_xyzo_min = m_transformInertial2This.transformReverse(bbox_xyz_min);
  m_bbox_xyzo_max = m_transformInertial2This.transformReverse(bbox_xyz_max);
}


void StructuredHexRaster::setupInterpolators()
{
  // nothing to be done here.
}


bool StructuredHexRaster::computeCCDataInterpolCoeffs(real x, real y, real z,
                                                 WeightedSet<real>& w_set)
{
  size_t ic_ref, jc_ref, kc_ref;
  bool inside;
  real w_octet[8];

  // May not request a value, if m_NumProtectLayers<1
  if(m_NumProtectLayers < 1) {
    cout << "StructuredHexRaster::computeCCDataInterpolCoeffs, m_NumProtectLayers = " << m_NumProtectLayers << endl;
  }

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
    w_set.clearWS();
    w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), w_octet[0]);
    w_set.pushBack(index(ic_ref,   jc_ref,   kc_ref+1), w_octet[1]);
    w_set.pushBack(index(ic_ref,   jc_ref+1, kc_ref  ), w_octet[2]);
    w_set.pushBack(index(ic_ref,   jc_ref+1, kc_ref+1), w_octet[3]);
    w_set.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ), w_octet[4]);
    w_set.pushBack(index(ic_ref+1, jc_ref,   kc_ref+1), w_octet[5]);
    w_set.pushBack(index(ic_ref+1, jc_ref+1, kc_ref  ), w_octet[6]);
    w_set.pushBack(index(ic_ref+1, jc_ref+1, kc_ref+1), w_octet[7]);
    //.. Eliminate round of contributors with below eps-weight
    /// @todo m_Eps is used twice, see also computeDeltas. Potentially conflicting.
    w_set.EliminateBelowEps(10*m_Eps, false, false);
    //.. Ensure sum of weights to be 1. to correct eps-errors
    w_set.adjustWeightSumShift(1.);
  }
  return inside;
}

bool StructuredHexRaster::computeCCGrad1NInterpolCoeffs(real x, real y, real z,
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

  // May not request a value, if m_NumProtectLayers<2
  if(m_NumProtectLayers<2) {
    cout << "StructuredHexRaster::computeCCGrad1NInterpolCoeffs, m_NumProtectLayers = " << m_NumProtectLayers << endl;
  }
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


void StructuredHexRaster::computeCCInterpolWeights(real& x, real& y, real& z,
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

  // StructuredHex distances to bounding element sides
  real x_ref = 0.5*m_Dx + ic_ref*m_Dx;
  real y_ref = 0.5*m_Dy + jc_ref*m_Dy;
  real z_ref = 0.5*m_Dz + kc_ref*m_Dz;

  computeInterpolWeights(x, y, z,
                         x_ref, y_ref, z_ref,
                         w_octet);
}

void StructuredHexRaster::computeNCInterpolWeights(real& x, real& y, real& z,
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

  // StructuredHex distances to bounding element sides
  real x_ref = in_ref*m_Dx;
  real y_ref = jn_ref*m_Dy;
  real z_ref = kn_ref*m_Dz;

  computeInterpolWeights(x, y, z,
                         x_ref, y_ref, z_ref,
                         w_octet);
}

void StructuredHexRaster::computeInterpolWeights(real& x, real& y, real& z,
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

  // StructuredHex distances to box sides
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

};


// void StructuredHexRaster::computeCCGrad1NInterpolCoeffs(real& x, real& y, real& z,
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

bool StructuredHexRaster::xyzToRefInterCell(real x, real y, real z,
				       size_t& ic_ref, size_t& jc_ref, size_t& kc_ref)
{
  bool inside = true;

  /// @todo Eps-stuff is quite ugly down here

  // check limits
  if(x < m_xCCInterMin-m_EpsDX) {ic_ref = 0; inside = false;};
  if(y < m_yCCInterMin-m_EpsDY) {jc_ref = 0; inside = false;};
  if(z < m_zCCInterMin-m_EpsDZ) {kc_ref = 0; inside = false;};
  if(x > m_xCCInterMax+m_EpsDX) {ic_ref = m_NumI-1; inside = false;};
  if(y > m_yCCInterMax+m_EpsDY) {jc_ref = m_NumJ-1; inside = false;};
  if(z > m_zCCInterMax+m_EpsDZ) {kc_ref = m_NumK-1; inside = false;};
  // get reference cell, if inside
  if(inside) {
    //.. Manipulate x,y,z if these coords depass the borders
    //   NOTE: inside==true means that x,y,z is in eps-bounds anyway
    //         if close, shift x,y,z to inside to ensure correct address pick
    if(x < m_xCCInterMin+m_EpsDX) {x = m_xCCInterMin+m_EpsDX;}
    if(y < m_yCCInterMin+m_EpsDY) {y = m_yCCInterMin+m_EpsDY;}
    if(z < m_zCCInterMin+m_EpsDZ) {z = m_zCCInterMin+m_EpsDZ;}
    if(x > m_xCCInterMax-m_EpsDX) {x = m_xCCInterMax-m_EpsDX;}
    if(y > m_yCCInterMax-m_EpsDY) {y = m_xCCInterMax-m_EpsDY;}
    if(z > m_zCCInterMax-m_EpsDZ) {z = m_xCCInterMax-m_EpsDZ;}
    //.. Position relative to lowest CC-coords.
    real ric_ref = (x - m_xCCMin) / m_Dx;
    real rjc_ref = (y - m_yCCMin) / m_Dy;
    real rkc_ref = (z - m_zCCMin) / m_Dz;
    //    iic_ref = lrint(floor(ric_ref));
    //    jjc_ref = lrint(floor(rjc_ref));
    //    kkc_ref = lrint(floor(rkc_ref));
    //.. Compute reference address. Ensure positive.
    //   NOTE: Due to Eps-errors, negative values for ric_ref, rjc_ref, rkc_ref are possible
    //         if m_NumProtectLayers==0 , probably never of interest, since no code will work
    //         without protection zone (?)
    ic_ref = size_t(ric_ref);
    jc_ref = size_t(rjc_ref);
    kc_ref = size_t(rkc_ref);
    //.. be sure to avoid size_t cycling. May never apply.
    if(ric_ref < 0.) {ic_ref = 0;}
    if(rjc_ref < 0.) {jc_ref = 0;}
    if(rkc_ref < 0.) {kc_ref = 0;}
    //.. Error checking, actually impossible: Be sure to be in addressing range
    bool error = false;
    if(ic_ref < m_NumProtectLayers) {
      cout << "StructuredHexRaster::xyzToRefInterCell, i_min" << endl; error = true;
    }
    if(jc_ref < m_NumProtectLayers) {
      cout << "StructuredHexRaster::xyzToRefInterCell, j_min" << endl; error = true;
    }
    if(kc_ref < m_NumProtectLayers) {
      cout << "StructuredHexRaster::xyzToRefInterCell, k_min" << endl; error = true;
    }
    if(ic_ref > m_NumI-2-m_NumProtectLayers) {
      cout << "StructuredHexRaster::xyzToRefInterCell, i_max" << endl; error = true;
    }
    if(jc_ref > m_NumJ-2-m_NumProtectLayers) {
      cout << "StructuredHexRaster::xyzToRefInterCell, j_max" << endl; error = true;
    }
    if(kc_ref > m_NumK-2-m_NumProtectLayers) {
      cout << "StructuredHexRaster::xyzToRefInterCell, k_max" << endl; error = true;
    }
    if(error) {
      /// @todo may put this in debug condition
      cout << endl;
      cout << " m_xyzCCInterMin = " << m_xCCInterMin << ", " << m_yCCInterMin << ", " << m_zCCInterMin << endl;
      cout << " m_xyzCCInterMax = " << m_xCCInterMax << ", " << m_yCCInterMax << ", " << m_zCCInterMax << endl;
      cout << " (x,y,z)         = " << x             << ", " << y             << ", " << z             << endl;
      BUG;
    }
  }
  return inside;
}

bool StructuredHexRaster::xyzToRefCell(real x, real y, real z,
				  size_t& ic_ref, size_t& jc_ref, size_t& kc_ref)
{
  bool inside = true;

  /// @todo Eps-stuff is quite ugly down here

  // check limits
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
    //.. Position relative to lowest CC-coords.
    real ric_ref = (x - m_xCCMin) / m_Dx;
    real rjc_ref = (y - m_yCCMin) / m_Dy;
    real rkc_ref = (z - m_zCCMin) / m_Dz;
    //    iic_ref = lrint(floor(ric_ref));
    //    jjc_ref = lrint(floor(rjc_ref));
    //    kkc_ref = lrint(floor(rkc_ref));
    //.. Compute reference address. Ensure positive.
    //   NOTE: Due to Eps-errors, negative values for ric_ref, rjc_ref, rkc_ref are possible
    //         if Eps_DXYZ==0.
    ic_ref = size_t(ric_ref);
    jc_ref = size_t(rjc_ref);
    kc_ref = size_t(rkc_ref);
    //.. Error checking, actually impossible: Be sure to be in addressing range
    bool error = false;
    if(ric_ref < 0.) {  // use real ric_ref, due to size_t cycling
      cout << "StructuredHexRaster::xyzToRefCell, i_min" << endl; error = true;
    }
    if(rjc_ref < 0.) {
      cout << "StructuredHexRaster::xyzToRefCell, j_min" << endl; error = true;
    }
    if(rkc_ref < 0.) {
      cout << "StructuredHexRaster::xyzToRefCell, k_min" << endl; error = true;
    }
    if(ic_ref > m_NumI-2) {
      cout << "StructuredHexRaster::xyzToRefCell, i_max" << endl; error = true;
    }
    if(jc_ref > m_NumJ-2) {
      cout << "StructuredHexRaster::xyzToRefCell, j_max" << endl; error = true;
    }
    if(kc_ref > m_NumK-2) {
      cout << "StructuredHexRaster::xyzToRefCell, k_max" << endl; error = true;
    }
    if(error) {
      /// @todo may put this in debug condition
      cout << endl;
      cout << " m_xyzCCMin = " << m_xCCMin << ", " << m_yCCMin << ", " << m_zCCMin << endl;
      cout << " m_xyzCCMax = " << m_xCCMax << ", " << m_yCCMax << ", " << m_zCCMax << endl;
      cout << " (x,y,z)    = " << x        << ", " << y        << ", " << z        << endl;
      BUG;
    }
  }
  return inside;
}

bool StructuredHexRaster::xyzToRefNode(real x, real y, real z,
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
    real rin_ref = (x - 0.) / m_Dx;
    real rjn_ref = (y - 0.) / m_Dy;
    real rkn_ref = (z - 0.) / m_Dz;
    //    iin_ref = lrint(floor(rin_ref));
    //    jjn_ref = lrint(floor(rjn_ref));
    //    kkn_ref = lrint(floor(rkn_ref));
    //.. Compute reference address. Ensure positive.
    //   NOTE: Due to Eps-errors, negative values for rin_ref, rjn_ref, rkn_ref are possible
    //         if Eps_DXYZ==0.
    in_ref = size_t(rin_ref);
    jn_ref = size_t(rjn_ref);
    kn_ref = size_t(rkn_ref);
    //.. Error checking, actually impossible: Be sure to be in addressing range
    bool error = false;
    if(rin_ref < 0.) {
      cout << "StructuredHexRaster::xyzToRefNode, i_min" << endl; error = true;
    }
    if(rjn_ref < 0.) {
      cout << "StructuredHexRaster::xyzToRefNode, j_min" << endl; error = true;
    }
    if(rkn_ref < 0.) {
      cout << "StructuredHexRaster::xyzToRefNode, k_min" << endl; error = true;
    }
    if(in_ref > m_NumI-1) {
      cout << "StructuredHexRaster::xyzToRefNode, i_max" << endl; error = true;
    }
    if(jn_ref > m_NumJ-1) {
      cout << "StructuredHexRaster::xyzToRefNode, j_max" << endl; error = true;
    }
    if(kn_ref > m_NumK-1) {
      cout << "StructuredHexRaster::xyzToRefNode, k_max" << endl; error = true;
    }
    if(error) {
      /// @todo may put this in debug condition
      cout << endl;
      cout << " xyz Min = " << 0. << ", " << 0. << ", " << 0. << endl;
      cout << " xyz Max = " << m_Lx << ", " << m_Ly << ", " << m_Lz << endl;
      cout << " (x,y,z)    = " << x        << ", " << y        << ", " << z             << endl;
      BUG;
    }
  }
  return inside;
}

