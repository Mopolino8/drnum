#include "cartesianpatch.h"

#ifdef WITH_VTK
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <QVector>
#endif

/** @todo proposed coord naming convention:
 *        (Xo,Yo,Zo), (xo,yo,zo) or vec3_t XYZo, xyzo     : coords in syst. of origin
 *        (X,Y,Z)   , (x,y,z)    or vec3_t XYZ , xyz      : coords in syst. of patch
 *        (XX,YY,ZZ), (xx,yy,zz) or vec3_t XXYYZZ, xxyyzz : coords in syst. of any foreign patch
 *
 *  ATTENTION: This differs essentially to actual naming, as aparently all coords are only
 *             in the XYZo-system. In proposed coord-syst any CartesianPatch will have its own
 *             aligned system, no matter how it is actually oriented.
 *             Define an origin, usefull at node (i,j,k) = (0,0,0)
 *             Define XYZo_shift as the position of (0,0,0) in coord-syst of origin
 */

CartesianPatch::CartesianPatch()
  : Patch()
{
  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
  //do allways with computeDeltas() m_Interpol_Initialized = false;
  /// @todo need a better eps-handling. Is there a minimum real?
  m_Eps = 1.e-6;
}

void CartesianPatch::computeDeltas()
{
  double Lx = sqrt(sqr(m_UX) + sqr(m_UY) + sqr(m_UZ));
  double Ly = sqrt(sqr(m_VX) + sqr(m_VY) + sqr(m_VZ));
  double Lz = sqrt(sqr(m_WX) + sqr(m_WY) + sqr(m_WZ));
  countFlops(15);
  countSqrts(3);

  m_DX = Lx/m_NumI;
  m_DY = Ly/m_NumJ;
  m_DZ = Lz/m_NumK;
  countFlops(3);

  m_InvDX = 1.0/m_DX;
  m_InvDY = 1.0/m_DX;
  m_InvDZ = 1.0/m_DX;
  countFlops(3);

  m_xCCMin = 0.5 * m_DX;
  m_yCCMin = 0.5 * m_DY;
  m_zCCMin = 0.5 * m_DZ;
  m_xCCMax = m_UX - 0.5 * m_DX;
  m_yCCMax = m_VY - 0.5 * m_DY;
  m_zCCMax = m_WZ - 0.5 * m_DZ;
  countFlops(6);

  m_xCCInterMin = (m_NumProtectLayers + 0.5) * m_DX;
  m_yCCInterMin = (m_NumProtectLayers + 0.5) * m_DY;
  m_zCCInterMin = (m_NumProtectLayers + 0.5) * m_DZ;
  m_xCCInterMax = m_UX - (m_NumProtectLayers + 0.5) * m_DX;
  m_yCCInterMax = m_VY - (m_NumProtectLayers + 0.5) * m_DY;
  m_zCCInterMax = m_WZ - (m_NumProtectLayers + 0.5) * m_DZ;
  countFlops(6);

  m_EpsDX = m_DX * m_Eps;
  m_EpsDY = m_DY * m_Eps;
  m_EpsDZ = m_DZ * m_Eps;
  countFlops(3);
}

void CartesianPatch::setupAligned(real x1, real y1, real z1, real x2, real y2, real z2)
{

  /** @todo Use a more general alignment method with trafo-matrix and internal coord system.
   *  I assume these are co-aligned here. */

  m_Xo = x1;
  m_Yo = y1;
  m_Zo = z1;

  m_UX = x2 - x1;
  m_UY = 0;
  m_UZ = 0;
  countFlops(1);

  m_VX = 0;
  m_VY = y2 - y1;
  m_VZ = 0;
  countFlops(1);

  m_WX = 0;
  m_WY = 0;
  m_WZ = z2 - z1;
  countFlops(1);

  computeDeltas();

  // position relative to inertial system
  vec3_t shift_inertial2this;
  shift_inertial2this[0] = -m_Xo;
  shift_inertial2this[1] = -m_Yo;
  shift_inertial2this[2] = -m_Zo;
  m_transformInertial2This.setVector(shift_inertial2this);

  /// @todo rotationm NOT implemented!!!

}

void CartesianPatch::resize(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  deleteData();
  Patch::resize(m_NumI*m_NumJ*m_NumK);
  computeDeltas();
}

void CartesianPatch::buildBoundingBox()
{
  vec3_t bbox_xyz_min(0., 0., 0.);
  vec3_t bbox_xyz_max(m_UX, m_VY, m_WZ);  /// @todo change names UX, VY, WZ
  m_bbox_xyzo_min = m_transformInertial2This.transformReverse(bbox_xyz_min);
  m_bbox_xyzo_max = m_transformInertial2This.transformReverse(bbox_xyz_max);
}

void CartesianPatch::extractReceiveCells()
{
  bool any_error = false;
  bool error;
  size_t cell_h;
  /** @todo check!! very prone to cut&paste errors.
    * Can do better than below looping: inserts many duplicates, that will be eliminated later.
    */
  for(size_t layer=0; layer<m_NumOverlapLayers; layer++) {
    // i_cell=[0, m_NumOverlapLayers-1] and i_cell=[m_NumI-1, m_NumI-1-m_NumOverlapLayers-1]
    for(size_t j_cell=0; j_cell<m_NumJ; j_cell++) {
      for(size_t k_cell=0; k_cell<m_NumK; k_cell++) {
        cell_h = save_index(layer, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        m_receive_cells.push_back(cell_h);
        cell_h = save_index(m_NumI-1-layer, j_cell, k_cell,
                            error);
        any_error = any_error || error;
        m_receive_cells.push_back(cell_h);
      }
    }
    // j_cell=[0, m_NumOverlapLayers-1] and j_cell=[m_NumJ-1, m_NumJ-1-m_NumOverlapLayers-1]
    for(size_t k_cell=0; k_cell<m_NumK; k_cell++) {
      for(size_t i_cell=0; i_cell<m_NumI; i_cell++) {
        cell_h = save_index(i_cell, layer, k_cell,
                            error);
        any_error = any_error || error;
        m_receive_cells.push_back(cell_h);
        cell_h = save_index(i_cell, m_NumJ-1-layer, k_cell,
                            error);
        any_error = any_error || error;
        m_receive_cells.push_back(cell_h);
      }
    }
    // k_cell=[0, m_NumOverlapLayers-1] and k_cell=[m_NumK-1, m_NumK-1-m_NumOverlapLayers-1]
    for(size_t i_cell=0; i_cell<m_NumI; i_cell++) {
      for(size_t j_cell=0; j_cell<m_NumJ; j_cell++) {
        cell_h = save_index(i_cell, j_cell, layer,
                            error);
        any_error = any_error || error;
        m_receive_cells.push_back(cell_h);
        cell_h = save_index(i_cell, j_cell, m_NumK-1-layer,
                            error);
        any_error = any_error || error;
        m_receive_cells.push_back(cell_h);
      }
    }
  }
}

void CartesianPatch::computeDependencies(const size_t& i_neighbour)
{
  Patch* neighbour_patch = m_neighbours[i_neighbour].first;
  CoordTransformVV trans = m_neighbours[i_neighbour].second;
  // patch bounding normal vectors in coord syst of donor
  vector<vec3_t> nnxyz_ijk;
  nnxyz_ijk.resize(3);
  nnxyz_ijk[0] = trans.transfree(vec3_t(1., 0., 0.));
  nnxyz_ijk[1] = trans.transfree(vec3_t(0., 1., 0.));
  nnxyz_ijk[2] = trans.transfree(vec3_t(0., 0., 1.));
  // Process all receiving cells and get interpolation sets
  //  - get own coeffs
  //  - transform into system of neighbour patch
  //  - interpolate there if receiving cell hits neighbours core region
  for(size_t ll_rc=0; ll_rc < m_receive_cells.size(); ll_rc++) {
    size_t l_rc = m_receive_cells[ll_rc];
    //    real x_rc = cellX(l_rc);
    //    real y_rc = cellY(l_rc);
    //    real z_rc = cellZ(l_rc);
    //    real xx_rc, yy_rc, zz_rc;
    vec3_t xyz_rc = xyzCell(l_rc);
    vec3_t xxyyzz_rc = trans.transform(xyz_rc);
    WeightedSet<real> w_set;
    if(m_InterpolateData) {
      if(neighbour_patch->computeCCDataInterpolCoeffs(xxyyzz_rc,
                                                      w_set)) {
        m_InterCoeffData_WS[i_neighbour]->push(ll_rc, l_rc, w_set);  // note: since l_rc are unique, no "add" is required
        m_receive_cell_data_hits[ll_rc]++;
      }
      if(m_InterpolateGrad1N) {
        int direction = boundingNVecDirection(l_rc); // note: trivial function, that ommits normal storage
        //        if(neighbour_patch->computeCCGrad1NInterpolCoeffs(xx_rc, yy_rc, zz_rc,
        //                                                          nnxyz_ijk[direction][0], nnxyz_ijk[direction][1], nnxyz_ijk[direction][2],
        //                                                          w_set)) {
        if(neighbour_patch->computeCCGrad1NInterpolCoeffs(xxyyzz_rc, nnxyz_ijk[direction],
                                                          w_set)) {
          m_InterCoeffGrad1N_WS[i_neighbour]->push(ll_rc, l_rc, w_set);
          m_receive_cell_grad1N_hits[ll_rc]++;
        }
      }
    }
  }
};


void CartesianPatch::setupInterpolators()
{
  // nothing to be done here.
  /** @note it is likely that setupInterpolators will be meaningful for other
      * patch types. If allways envoqued, shift stuff from computeDeltas back here.
      */
}

// void CartesianPatch::setupInterpolators(size_t num_protectlayers)
// {

//   m_NumProtectLayers = num_protectlayers;

//   // Cell center coords of lowest address cell (i,j,k) = (0,0,0)
//   //  m_xCCMin = 0.5 * m_DX;
//   //  m_yCCMin = 0.5 * m_DY;
//   //  m_zCCMin = 0.5 * m_DZ;

//   // CC-coordinates limits to access data from foreign patches
//   m_xCCInterMin = (m_NumProtectLayers + 0.5) * m_DX;
//   m_yCCInterMin = (m_NumProtectLayers + 0.5) * m_DY;
//   m_zCCInterMin = (m_NumProtectLayers + 0.5) * m_DZ;
//   m_xCCInterMax = m_UX - (m_NumProtectLayers + 0.5) * m_DX;
//   m_yCCInterMax = m_VY - (m_NumProtectLayers + 0.5) * m_DY;
//   m_zCCInterMax = m_WZ - (m_NumProtectLayers + 0.5) * m_DZ;

//   m_EpsDX = m_Eps * m_DX;
//   m_EpsDY = m_Eps * m_DY;
//   m_EpsDZ = m_Eps * m_DZ;

//   //  m_Interpol_Initialized = true;
// }

bool CartesianPatch::computeCCDataInterpolCoeffs(real x, real y, real z,
                                                 WeightedSet<real>& w_set)
{
  size_t ic_ref, jc_ref, kc_ref;
  bool inside;
  real w_octet[8];

  // May not request a value, if m_NumProtectLayers<1
  if(m_NumProtectLayers < 1) {
    cout << "CartesianPatch::computeCCDataInterpolCoeffs, m_NumProtectLayers = " << m_NumProtectLayers << endl;
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

  // May not request a value, if m_NumProtectLayers<2
  if(m_NumProtectLayers<2) {
    cout << "CartesianPatch::computeCCGrad1NInterpolCoeffs, m_NumProtectLayers = " << m_NumProtectLayers << endl;
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
    real div_x = nnx/(2*m_DX);
    real div_y = nny/(2*m_DY);
    real div_z = nnz/(2*m_DZ);

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

  // cartesian distances to bounding element sides
  real x_ref = 0.5*m_DX + ic_ref*m_DX;
  real y_ref = 0.5*m_DY + jc_ref*m_DY;
  real z_ref = 0.5*m_DZ + kc_ref*m_DZ;

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

  // cartesian distances to bounding element sides
  real x_ref = in_ref*m_DX;
  real y_ref = jn_ref*m_DY;
  real z_ref = kn_ref*m_DZ;

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
  //   1:          ( x_ref,      y_ref,      z_ref+m_DZ )
  //   2:          ( x_ref,      y_ref+m_DY, z_ref      )
  //   3:          ( x_ref,      y_ref+m_DY, z_ref+m_DZ )
  //   4:          ( x_ref+m_DX, y_ref,      z_ref      )
  //   5:          ( x_ref+m_DX, y_ref,      z_ref+m_DZ )
  //   6:          ( x_ref+m_DX, y_ref+m_DY, z_ref      )
  //   7:          ( x_ref+m_DX, y_ref+m_DY, z_ref+m_DZ )

  // cartesian distances to box sides
  real low_diff_x_ref = x - x_ref;
  real low_diff_y_ref = y - y_ref;
  real low_diff_z_ref = z - z_ref;
  real up_diff_x_ref = m_DX - low_diff_x_ref;
  real up_diff_y_ref = m_DY - low_diff_y_ref;
  real up_diff_z_ref = m_DZ - low_diff_z_ref;

  real inv_cell_vol = m_InvDX * m_InvDY * m_InvDZ;

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
// d_dx_04.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDX);
// d_dx_04.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDX);
// d_dx_04 *= up_diff_y_ref;
// d_dx_04 *= up_diff_z_ref;
// //.. cell-cell difference grad for 1->5
// WeightedSet<real> d_dx_15();
// d_dx_15.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDX);
// d_dx_15.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDX);
// d_dx_15 *= up_diff_y_ref;
// d_dx_15 *= low_diff_z_ref;
// //.. cell-cell difference grad for 2->6
// WeightedSet<real> d_dx_26();
// d_dx_26.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDX);
// d_dx_26.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDX);
// d_dx_26 *= low_diff_y_ref;
// d_dx_26 *= up_diff_z_ref;
// //.. cell-cell difference grad for 3->7
// WeightedSet<real> d_dx_37();
// d_dx_37.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDX);
// d_dx_37.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDX);
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
// d_dx_02.pushBack(index(ic_ref,   jc_ref,   kc_ref  ), -m_InvDY);
// d_dx_02.pushBack(index(ic_ref+1, jc_ref,   kc_ref  ),  m_InvDY);
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
// face_i_low *= (m_InvDY * m_InvDZ);
// face_i_up  *= (m_InvDY * m_InvDZ);
// WeightedSet<real> d_dx();
// d_dx.clearWS();
// d_dx.Concatenate(face_i_up ,  m_InvDX);
// d_dx.Concatenate(face_i_low, -m_InvDX);

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
// face_j_low *= (m_InvDX * m_InvDZ);
// face_j_up  *= (m_InvDX * m_InvDZ);
// WeightedSet<real> d_dx();
// d_dx.clearWS();
// d_dx.Concatenate(face_j_up ,  m_InvDY);
// d_dx.Concatenate(face_j_low, -m_InvDY);

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
// face_i_low *= (m_InvDY * m_InvDZ);
// face_i_up  *= (m_InvDY * m_InvDZ);
// WeightedSet<real> d_dx();
// d_dx.clearWS();
// d_dx.Concatenate(face_i_up ,  m_InvDX);
// d_dx.Concatenate(face_i_low, -m_InvDX);

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

bool CartesianPatch::xyzToRefInterCell(real x, real y, real z,
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
    real ric_ref = (x - m_xCCMin) / m_DX;
    real rjc_ref = (y - m_yCCMin) / m_DY;
    real rkc_ref = (z - m_zCCMin) / m_DZ;
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
    if(ic_ref < m_NumProtectLayers) {cout << "CartesianPatch::xyzToRefInterCell, i_min"; error = true;}
    if(jc_ref < m_NumProtectLayers) {cout << "CartesianPatch::xyzToRefInterCell, j_min"; error = true;}
    if(kc_ref < m_NumProtectLayers) {cout << "CartesianPatch::xyzToRefInterCell, k_min"; error = true;}
    if(ic_ref > m_NumI-2-m_NumProtectLayers) {cout << "CartesianPatch::xyzToRefInterCell, i_max"; error = true;}
    if(jc_ref > m_NumJ-2-m_NumProtectLayers) {cout << "CartesianPatch::xyzToRefInterCell, j_max"; error = true;}
    if(kc_ref > m_NumK-2-m_NumProtectLayers) {cout << "CartesianPatch::xyzToRefInterCell, k_max"; error = true;}
    if(error) {
      /// @todo may put this in debug condition
      cout << endl;
      cout << " m_xyzCCInterMin = " << m_xCCInterMin << ", " << m_yCCInterMin << ", " << m_zCCInterMin << endl;
      cout << " m_xyzCCInterMax = " << m_xCCInterMax << ", " << m_yCCInterMax << ", " << m_zCCInterMax << endl;
      cout << " (x,y,z)         = " << x             << ", " << y             << ", " << z             << endl;
      exit(EXIT_FAILURE);
    }
  }
  return inside;
}

bool CartesianPatch::xyzToRefCell(real x, real y, real z,
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
    if(y > m_yCCMax-m_EpsDY) {y = m_xCCMax-m_EpsDY;}
    if(z > m_zCCMax-m_EpsDZ) {z = m_xCCMax-m_EpsDZ;}
    //.. Position relative to lowest CC-coords.
    real ric_ref = (x - m_xCCMin) / m_DX;
    real rjc_ref = (y - m_yCCMin) / m_DY;
    real rkc_ref = (z - m_zCCMin) / m_DZ;
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
      cout << "CartesianPatch::xyzToRefInterCell, i_min"; error = true;
      ic_ref = 0;
    }
    if(rjc_ref < 0.) {
      cout << "CartesianPatch::xyzToRefInterCell, j_min"; error = true;
      jc_ref = 0;
    }
    if(rkc_ref < 0.) {
      cout << "CartesianPatch::xyzToRefInterCell, k_min"; error = true;
      kc_ref = 0;
    }
    if(ic_ref > m_NumI-2) {
      cout << "CartesianPatch::xyzToRefInterCell, i_max"; error = true;
    }
    if(jc_ref > m_NumJ-2) {
      cout << "CartesianPatch::xyzToRefInterCell, j_max"; error = true;
    }
    if(kc_ref > m_NumK-2) {
      cout << "CartesianPatch::xyzToRefInterCell, k_max"; error = true;
    }
    if(error) {
      /// @todo may put this in debug condition
      cout << endl;
      cout << " m_xyzCCMin = " << m_xCCMin << ", " << m_yCCMin << ", " << m_zCCMin << endl;
      cout << " m_xyzCCMax = " << m_xCCMax << ", " << m_yCCMax << ", " << m_zCCMax << endl;
      cout << " (x,y,z)    = " << x        << ", " << y        << ", " << z        << endl;
      exit(EXIT_FAILURE);
    }
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
  if(x > m_UX+m_EpsDX) {in_ref = m_NumI-1; inside = false;};
  if(y > m_VY+m_EpsDY) {jn_ref = m_NumJ-1; inside = false;};
  if(z > m_WZ+m_EpsDZ) {kn_ref = m_NumK-1; inside = false;};

  // get reference cell, if inside
  if(inside) {
    //.. Manipulate x,y,z if these coords depass the borders
    //   NOTE: inside==true means that x,y,z is in eps-bounds anyway
    //         if close, shift x,y,z to inside to ensure correct address pick
    if(x < m_xCCMin+m_EpsDX) {x = m_xCCMin+m_EpsDX;}
    if(y < m_yCCMin+m_EpsDY) {y = m_yCCMin+m_EpsDY;}
    if(z < m_zCCMin+m_EpsDZ) {z = m_zCCMin+m_EpsDZ;}
    if(x > m_xCCMax-m_EpsDX) {x = m_xCCMax-m_EpsDX;}
    if(y > m_yCCMax-m_EpsDY) {y = m_xCCMax-m_EpsDY;}
    if(z > m_zCCMax-m_EpsDZ) {z = m_xCCMax-m_EpsDZ;}
    //.. Position relative to lowest CC-coords.
    real rin_ref = (x - m_xCCMin) / m_DX;
    real rjn_ref = (y - m_yCCMin) / m_DY;
    real rkn_ref = (z - m_zCCMin) / m_DZ;
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
      cout << "CartesianPatch::xyzToRefInterCell, i_min"; error = true;
      in_ref = 0;
    }
    if(rjn_ref < 0.) {
      cout << "CartesianPatch::xyzToRefInterCell, j_min"; error = true;
      jn_ref = 0;
    }
    if(rkn_ref < 0.) {
      cout << "CartesianPatch::xyzToRefInterCell, k_min"; error = true;
      kn_ref = 0;
    }
    if(in_ref > m_NumI-2) {
      cout << "CartesianPatch::xyzToRefInterCell, i_max"; error = true;
    }
    if(jn_ref > m_NumJ-2) {
      cout << "CartesianPatch::xyzToRefInterCell, j_max"; error = true;
    }
    if(kn_ref > m_NumK-2) {
      cout << "CartesianPatch::xyzToRefInterCell, k_max"; error = true;
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


#ifdef WITH_VTK

void CartesianPatch::writeToVtk(QString file_name)
{
  vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();

  vtkSmartPointer<vtkFloatArray> xc = vtkSmartPointer<vtkFloatArray>::New();
  for (size_t i = 0; i < m_NumI + 1; ++i) {
    xc->InsertNextValue(i*dx());
  }

  vtkSmartPointer<vtkFloatArray> yc = vtkSmartPointer<vtkFloatArray>::New();
  for (size_t j = 0; j < m_NumJ + 1; ++j) {
    yc->InsertNextValue(j*dy());
  }

  vtkSmartPointer<vtkFloatArray> zc = vtkSmartPointer<vtkFloatArray>::New();
  for (size_t k = 0; k < m_NumK + 1; ++k) {
    zc->InsertNextValue(k*dz());
  }

  grid->SetDimensions(m_NumI + 1, m_NumJ + 1, m_NumK + 1);
  grid->SetXCoordinates(xc);
  grid->SetYCoordinates(yc);
  grid->SetZCoordinates(zc);

  for (size_t i_field = 0; i_field < numFields(); ++i_field) {
    QString field_name;
    field_name.setNum(i_field + 1);
    if (field_name.size() < 2) {
      field_name = "0" + field_name;
    }
    for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
      QString var_name;
      var_name.setNum(i_var + 1);
      if (var_name.size() < 2) {
        var_name = "0" + var_name;
      }
      vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
      var->SetName(qPrintable("f_" + field_name + "_" + var_name));
      var->SetNumberOfValues(variableSize());
      grid->GetCellData()->AddArray(var);
      vtkIdType id = 0;
      for (size_t k = 0; k < m_NumK; ++k) {
        for (size_t j = 0; j < m_NumJ; ++j) {
          for (size_t i = 0; i < m_NumI; ++i) {
            var->SetValue(id, f(getVariable(i_field, i_var), i, j, k));
            ++id;
          }
        }
      }
    }
  }

  vtkSmartPointer<vtkXMLRectilinearGridWriter> vtr = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
  vtr->SetFileName(qPrintable(file_name + ".vtr"));
  vtr->SetInput(grid);
  vtr->Write();
}

#endif
