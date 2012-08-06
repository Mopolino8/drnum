#ifndef CARTESIANPATCH_H
#define CARTESIANPATCH_H

class CartesianPatch;

#include "patch.h"

#ifdef WITH_VTK
#include <QString>
#endif


class CartesianPatch : public Patch
{

protected: // attributes

  size_t m_NumI;
  size_t m_NumJ;
  size_t m_NumK;

  /// @todo Propose all cartesian patches with internal coord-syst with origin at node (0,0,0)

  real m_Xo;   /// @todo propose this to be allways 0.
  real m_Yo;   /// @todo propose this to be allways 0.
  real m_Zo;   /// @todo propose this to be allways 0.
  real m_UX;
  real m_UY;
  real m_UZ;
  real m_VX;
  real m_VY;
  real m_VZ;
  real m_WX;
  real m_WY;
  real m_WZ;
  real m_DX;
  real m_DY;
  real m_DZ;
  real m_InvDX;
  real m_InvDY;
  real m_InvDZ;

  real m_bbox_xyz_min;
  real m_bbox_xyz_max;
  bool m_bbox_OK;

  real m_LimiterEpsilon;

  /// @todo change to vec3_t later
  //  SVec3 m_xyzCCMin;      ///< Cell center coords of lowest address cell
  //  SVec3 m_xyzInterMin;   ///< lower CC coord-limits to access data from foreign patches
  //  SVec3 m_xyzInterMax;   ///< upper CC coord-limits to access data from foreign patches
  //->patch  size_t m_NumProtectLayers;    ///< Number of protection layers
  real m_xCCMin,      m_yCCMin,      m_zCCMin;
  real m_xCCMax,      m_yCCMax,      m_zCCMax;
  real m_xCCInterMin, m_yCCInterMin, m_zCCInterMin;
  real m_xCCInterMax, m_yCCInterMax, m_zCCInterMax;
  real m_Eps, m_EpsDX, m_EpsDY, m_EpsDZ;

protected: // methods

  void computeDeltas();
  virtual void extractReceiveCells();
  virtual void buildBoundingBox();

public: // methods

  CartesianPatch();
  void setupAligned(real x1, real y1, real z1, real x2, real y2, real z2);
  void resize(size_t num_i, size_t num_j, size_t num_k);
  virtual void computeDependencies(const size_t& i_neighbour);

  size_t sizeI() { return m_NumI; }
  size_t sizeJ() { return m_NumJ; }
  size_t sizeK() { return m_NumK; }


  /// @todo shift all addressing stuff into StructuredPatch and inherite CartesianPatch from it
  /// @todo I'm not sure, if the int - size_t conversion might cause performance issues

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the index in the one dimensional data field
   */
  size_t index(int i, int j, int k) { return i*m_NumJ*m_NumK + j*m_NumK + k; }

  /**
   * @brief Get the indicees (i, j, k) from field index/
   * @param l the field index
   * @param i first Cartesian index  (return reference)
   * @param j second Cartesian index (return reference)
   * @param k third Cartesian index  (return reference)
   * @return the index in the one dimensional data field
   */
  inline void ijk(const size_t& l,
                  size_t& i, size_t& j, size_t& k) const {
    //    div_t divresult;
    //    divresult = div (l, m_NumJ*m_NumK);  // causes compile error (ambiguous)
    //    i = divresult.quot;
    //    divresult = div (divresult.rem, m_NumK);
    //    j = divresult.quot;
    //    k = divresult.rem;
    size_t rest;
    i = l / (m_NumJ*m_NumK);
    rest = l - i*(m_NumJ*m_NumK);
    j = rest / m_NumK;
    k = rest - j*m_NumK;
  }

  inline vec3_t xyzCell(const size_t& l) const
  {
    size_t i, j, k;
    ijk(l,
        i, j, k);
    return xyzCell(i, j, k);
  }

  inline vec3_t xyzCell(const size_t& i, const size_t& j, const size_t& k) const
  {
    vec3_t xyz;
    xyzCell(i, j, k,
            xyz[0], xyz[1], xyz[2]);
    return xyz;
  }

  inline void xyzCell(const size_t& i, const size_t& j, const size_t& k,
                      real& x, real& y, real& z) const
  {
    x = m_xCCMin + i*m_DX;
    y = m_yCCMin + j*m_DY;
    z = m_zCCMin + k*m_DZ;
  }

  inline int boundingNVecDirection(const size_t& lc) const
  {
    int direction;
    size_t i, j, k;
    size_t min_index_dist, min_index_dist_h;
    ijk(lc,
        i, j, k);
    // prime with i-direction: 0
    min_index_dist = min(i, (m_NumI-1-i));
    direction = 0;
    // check j-direction: 1
    min_index_dist_h = min(j, (m_NumJ-1-j));
    if(min_index_dist > min_index_dist_h) {
      min_index_dist = min_index_dist_h;
      direction = 1;
    }
    min_index_dist_h = min(k, (m_NumK-1-k));
    if(min_index_dist > min_index_dist_h) {
      min_index_dist = min_index_dist_h;
      direction = 2;
    }
    return direction;
  }

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param error true, if (i, j, k) are out of bounds
   * @return the index in the one dimensional data field
   */
  size_t save_index(int i, int j, int k,
                    bool& error) {
    size_t si = i;  // avoid vast compiler warnings
    size_t sj = j;
    size_t sk = k;
    error = false;
    if(i < 0) error = true;
    if(j < 0) error = true;
    if(k < 0) error = true;
    if(si > m_NumI-1) error = true;
    if(sj > m_NumJ-1) error = true;
    if(sk > m_NumK-1) error = true;
    return si*m_NumJ*m_NumK + sj*m_NumK + sk;
  }

  /**
   * @brief Get the value of a variable at an (i, j, k) triple.
   * @param field a pointer to the variable data
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the field value at (i, j, k).
   */
  real& f(real *var, size_t i, size_t j, size_t k) { return var[i*m_NumJ*m_NumK + j*m_NumK + k]; }

  /**
   * @brief Get the value of a variable at an (i, j, k) triple.
   * @param i_field field index
   * @param i_var variable index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the field value at (i, j, k).
   */
  real& f(size_t i_field, size_t i_var, size_t i, size_t j, size_t k) { return getVariable(i_field, i_var)[i*m_NumJ*m_NumK + j*m_NumK + k]; }

  real dx() { return m_DX; }
  real dy() { return m_DY; }
  real dz() { return m_DZ; }
  real dV() { return m_DX*m_DY*m_DZ; }
  real idx() { return m_InvDX; }
  real idy() { return m_InvDY; }
  real idz() { return m_InvDZ; }


  /**
   * Set up interpolation methods for giving data to foreign patches.
   * Example: Build up Split- or Octrees for search operations, etc ... depending on patch type.
   * @param num_protection number of overlap cell layers in which no data access is permissible.
   */
  virtual void setupInterpolators();

  /**
   * Get data interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCDataInterpolCoeffs(real x, real y, real z,
                                           WeightedSet<real>& w_set);

  /**
   * Get directional derivative (grad*n) interpolation coeff-sets.
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param nx the x-component of directional vector in the coords of the present patch
   * @param ny the y-component of directional vector in the coords of the present patch
   * @param nz the z-component of directional vector in the coords of the present patch
   * @param w_set WeightedSet<real> object on which to write data (return reference)
   * @return true, if interpol sets were found.
   */
  virtual bool computeCCGrad1NInterpolCoeffs(real x, real y, real z,
					     real nx, real ny, real nz,
                                             WeightedSet<real>& w_set);

  /**
   * Compute a set of CC interpolation weights for a given point (x,y,z) in inverse mesh.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param ic_ref i-index of reference cell
   * @param jc_ref j-index of reference cell
   * @param kc_ref k-index of reference cell
   * @param w_octet eight values, that represent interpolation weights for cells in intervall (return reference)
   */
  void computeCCInterpolWeights(real& x, real& y, real& z,
                                size_t& ic_ref, size_t& jc_ref, size_t& kc_ref,
                                real* w_octet);

  /**
   * Compute a set of NC interpolation weights for a given point (x,y,z) in natural mesh.
   * NOTE: mesh cell formed by nodes in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param in_ref i-index of reference node
   * @param jn_ref j-index of reference node
   * @param kn_ref k-index of reference node
   * @param w_octet eight values, that represent interpolation weights for cells in intervall (return reference)
   */
  void computeNCInterpolWeights(real& x, real& y, real& z,
                                size_t& in_ref, size_t& jn_ref, size_t& kn_ref,
                                real* w_octet);

  /**
   * Compute a set of interpolation weights for a given point (x,y,z) in a box.
   * NOTE: box corner coords in 3D intervall
   * ([x_ref..x_ref+m_DX], [y_ref..y_ref+m_DY], [z_ref..z_ref+m_DZ]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param in_ref i-index of reference node
   * @param jn_ref j-index of reference node
   * @param kn_ref k-index of reference node
   * @param w_octet eight values, that represent interpolation weights for cells in intervall (return reference)
   */
  void computeInterpolWeights(real& x, real& y, real& z,
                              real& x_ref, real& y_ref, real& z_ref,
                              real* w_octet);

  /**
   * Find inverse mesh cell, in which the point (x,y,z) resides.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param i_ref i-index of reference cell (return reference)
   * @param j_ref j-index of reference cell (return reference)
   * @param k_ref k-index of reference cell (return reference)
   * @return boolean, indicating point (x,y,z) is inside core region
   */
  bool xyzToRefInterCell(real x, real y, real z,
			 size_t& i_ref, size_t& j_ref, size_t& k_ref);

  /**
   * Find inverse mesh cell, in which the point (x,y,z) resides.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param ic_ref i-index of reference cell (return reference)
   * @param jc_ref j-index of reference cell (return reference)
   * @param kc_ref k-index of reference cell (return reference)
   * @return boolean, indicating point (x,y,z) is inside core region
   */
  bool xyzToRefCell(real x, real y, real z,
                    size_t& ic_ref, size_t& jc_ref, size_t& kc_ref);

  /**
   * Find cell, in which the point (x,y,z) resides.
   * NOTE: inverse mesh cell formed by cell centers in 3D address-intervall
   * ([i_ref..i_ref+1], [j_ref..j_ref+1], [k_ref..k_ref+1]).
   * @param x the x-value in the coords of the present patch
   * @param y the y-value in the coords of the present patch
   * @param z the z-value in the coords of the present patch
   * @param in_ref i-index of reference node (return reference)
   * @param jn_ref j-index of reference node (return reference)
   * @param kn_ref k-index of reference node (return reference)
   * @return boolean, indicating point (x,y,z) is inside core region
   */
  bool xyzToRefNode(real x, real y, real z,
                    size_t& in_ref, size_t& jn_ref, size_t& kn_ref);

#ifdef WITH_VTK
  void writeToVtk(QString file_name);
#endif

};


#endif // CARTESIANPATCH_H
