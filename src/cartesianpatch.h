#ifndef CARTESIANPATCH_H
#define CARTESIANPATCH_H

class CartesianPatch;

#include "patch.h"

#ifdef WITH_VTK
#include <QString>
#include <QVector>
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#endif

class CartesianPatch : public Patch
{

protected: // attributes

  size_t m_NumI; // number of CELLS in I-direction (num nodes = m_NumI+1)
  size_t m_NumJ; // number of CELLS in J-direction (num nodes = m_NumJ+1)
  size_t m_NumK; // number of CELLS in K-direction (num nodes = m_NumK+1)

//  real m_Xo; // put reference position in patch.h
//  real m_Yo; // put reference position in patch.h
//  real m_Zo; // put reference position in patch.h

  /** @todo ommit the following, if no longer needed (setupAligned) */
  real m_Uxo;
  real m_Uyo;
  real m_Uzo;

  real m_Vxo;
  real m_Vyo;
  real m_Vzo;

  real m_Wxo;
  real m_Wyo;
  real m_Wzo;
  /// until here

  real m_Lx;
  real m_Ly;
  real m_Lz;

  real m_Dx;
  real m_Dy;
  real m_Dz;
  real m_InvDx;
  real m_InvDy;
  real m_InvDz;

  size_t m_numProtXmin;
  size_t m_numProtXmax;
  size_t m_numProtYmin;
  size_t m_numProtYmax;
  size_t m_numProtZmin;
  size_t m_numProtZmax;

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

  /**
   * Constructor
   * NOTE: The default number of protection as well as overplap layers is 1
   * @param num_protectlayers number of protection layers, in which no foreign data access is allowed
   * @param num_overlaplayers number of overlap layers, for which to get data from donor neighbour patches
   */
  CartesianPatch(size_t num_protectlayers=1, size_t num_overlaplayers=1);

  /**
    * Set number of protection layers on all boundaries of the CartesianPatch.
    * NOTE: method setNumProtectException(...) allows individual settings for all six boundaries of CartesianPatch.
    * @param num_protectlayers number of protection layers
    */
  virtual void setNumProtectLayers(size_t num_protectlayers);

  /**
    * Read mesh data from file
    * @param s_mesh the stream to read from
    * @return true, if successful
    */
  virtual bool readFromFile(ifstream& s_mesh);

  /**
    * Write mesh data to file
    * @param s_mesh the stream to write to
    * @return true, if successful
    */
  virtual bool writeToFile(ifstream& s_mesh);

  /**
    * Scale patch relative to origin of parental coordsyst.
    * NOTE: Affects reference position and physical patch size.
    * Virtual: base class Patch holds reference position. Derived class holds phys. patch size.
    * @param scfactor scaling factor.
    */
  virtual void scaleRefParental(real scfactor);

  /**
    * Set individual number of protection layers on all six boundaries of the CartesianPatch. This setting is
    * intended to allow exceptions, like interpolations along boundaries or reduced (1D, 2D) computations.
    * @param numProtXmin number of protection layers on Xmin-side
    * @param numProtXmax number of protection layers on Xmax-side
    * @param numProtYmin number of protection layers on Ymin-side
    * @param numProtYmax number of protection layers on Ymax-side
    * @param numProtZmin number of protection layers on Zmin-side
    * @param numProtZmax number of protection layers on Zmax-side
    */
  void setNumProtectException(const size_t& numProtXmin, const size_t& numProtXmax,
                              const size_t& numProtYmin, const size_t& numProtYmax,
                              const size_t& numProtZmin, const size_t& numProtZmax);

  /** @TODO: Eliminate this method and use setupMetrics instead, if working correctly. */
  void setupAligned(real xo1, real yo1, real zo1, real xo2, real yo2, real zo2);

  /**
    * Set up the metrics of the block.
    * @param ilength physical size of block in i direction
    * @param jlength physical size of block in j direction
    * @param klength physical size of block in k direction
    */
  void setupMetrics(real ilength, real jlength, real klength);
//  void setupMetrics(real xo1, real yo1, real zo1,
//                    real ilength, real jlength, real klength,
//                    real base_ix, base_iy, base_ik,
//                    real base_jx, base_jy, base_jk);

  void resize(size_t num_i, size_t num_j, size_t num_k);
  virtual bool computeDependencies(const size_t& i_neighbour);

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
    x = m_xCCMin + i*m_Dx;
    y = m_yCCMin + j*m_Dy;
    z = m_zCCMin + k*m_Dz;
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
   * @brief Get a variable set at a specified (i,j,k) position.
   * @param i_field the field index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param var will hold the conservative variable set afterwards (needs to be allocated beforehand)
   */
  void getVar(size_t i_field, size_t i, size_t j, size_t k, real* var);

  /**
   * @brief Set a variable set at a specified (i,j,k) position.
   * @param i_field the field index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param the conservative variable set
   */
  void setVar(size_t i_field, size_t i, size_t j, size_t k, real* var);

  /**
   * @brief Get the value of a variable at an (i, j, k) triple.
   * @param i_field field index
   * @param i_var variable index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the field value at (i, j, k).
   */
  real& f(size_t i_field, size_t i_var, size_t i, size_t j, size_t k);

  /**
   * Get the gradient in x direction at a specifed (i,j,k) position.
   * This method will automaitically respect domain borders (one sided gradient)
   * @param i_field the field index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param grad will hold the conservative x gradient set afterwards (needs to be allocated beforehand)
   */
  void getXGrad(size_t i_field, size_t i, size_t j, size_t k, real* grad);

  /**
   * Get the gradient in y direction at a specifed (i,j,k) position.
   * This method will automaitically respect domain borders (one sided gradient)
   * @param i_field the field index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param grad will hold the conservative y gradient set afterwards (needs to be allocated beforehand)
   */
  void getYGrad(size_t i_field, size_t i, size_t j, size_t k, real* grad);

  /**
   * Get the gradient in z direction at a specifed (i,j,k) position.
   * This method will automaitically respect domain borders (one sided gradient)
   * @param i_field the field index
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @param grad will hold the conservative z gradient set afterwards (needs to be allocated beforehand)
   */
  void getZGrad(size_t i_field, size_t i, size_t j, size_t k, real* grad);

  /**
   * Check if an (i,j,k) triple is inside the patch.
   * Attention only the upper limit will be checked (unsigned data type).
   * @param i first index
   * @param j second index
   * @param k third index
   * @return true if it is a valid (i,j,k) triple, false otherwise
   */
  bool checkRange(size_t i, size_t j, size_t k);

  real dx() { return m_Dx; }
  real dy() { return m_Dy; }
  real dz() { return m_Dz; }
  real dV() { return m_Dx*m_Dy*m_Dz; }
  real idx() { return m_InvDx; }
  real idy() { return m_InvDy; }
  real idz() { return m_InvDz; }

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
   * ([x_ref..x_ref+m_Dx], [y_ref..y_ref+m_Dy], [z_ref..z_ref+m_Dz]).
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
  template <typename TVariables> void writeToVtk(size_t i_field, QString file_name, const TVariables& variables);
#endif

};

inline real& CartesianPatch::f(size_t i_field, size_t i_var, size_t i, size_t j, size_t k)
{
#ifdef DEBUG
  if (!checkRange(i, j, k)) {
    BUG;
  }
#endif
  GlobalDebug::ijk(i, j, k);
  return getVariable(i_field, i_var)[i*m_NumJ*m_NumK + j*m_NumK + k];
}

inline void CartesianPatch::getVar(size_t i_field, size_t i, size_t j, size_t k, real *var)
{
  GlobalDebug::ijk(i, j, k);
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    var[i_var] = f(i_field, i_var, i, j, k);
  }
}

inline void CartesianPatch::setVar(size_t i_field, size_t i, size_t j, size_t k, real *var)
{
  GlobalDebug::ijk(i, j, k);
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    f(i_field, i_var, i, j, k) = var[i_var];
  }
}

inline void CartesianPatch::getXGrad(size_t i_field, size_t i, size_t j, size_t k, real *grad)
{
  real* var = new real[numVariables()];
  real D = 1.0/dx();
  countFlops(1);
  if (i > 0 && i < sizeI()-1) {
    getVar(i_field, i+1, j, k, grad);
    getVar(i_field, i-1, j, k, var);
    D *= 0.5;
    countFlops(1);
  } else if (i > 0) {
    getVar(i_field, i, j, k, grad);
    getVar(i_field, i-1, j, k, var);
  } else {
    getVar(i_field, i+1, j, k, grad);
    getVar(i_field, i, j, k, var);
  }
  for (int i_var = 0; i_var < numVariables(); ++i_var) {
    grad[i_var] -= var[i_var];
    grad[i_var] *= D;
  }
  countFlops(2*numVariables());
  delete [] var;
}

inline void CartesianPatch::getYGrad(size_t i_field, size_t i, size_t j, size_t k, real *grad)
{
  if (sizeJ() > 2) {
    real* var = new real[numVariables()];
    real D = 1.0/dy();
    countFlops(1);
    if (j > 0 && j < sizeJ()-1) {
      getVar(i_field, i, j+1, k, grad);
      getVar(i_field, i, j-1, k, var);
      D *= 0.5;
      countFlops(1);
    } else if (j > 0) {
      getVar(i_field, i, j, k, grad);
      getVar(i_field, i, j-1, k, var);
    } else {
      getVar(i_field, i, j+1, k, grad);
      getVar(i_field, i, j, k, var);
    }
    for (int i_var = 0; i_var < numVariables(); ++i_var) {
      grad[i_var] -= var[i_var];
      grad[i_var] *= D;
    }
    countFlops(2*numVariables());
    delete [] var;
  } else {
    for (int i_var = 0; i_var < numVariables(); ++i_var) {
      grad[i_var] = 0;
    }
  }
}

inline void CartesianPatch::getZGrad(size_t i_field, size_t i, size_t j, size_t k, real *grad)
{
  real* var = new real[numVariables()];
  real D = 1.0/dz();
  countFlops(1);
  if (k > 0 && k < sizeK()-1) {
    getVar(i_field, i, j, k+1, grad);
    getVar(i_field, i, j, k-1, var);
    D *= 0.5;
    countFlops(1);
  } else if (k > 0) {
    getVar(i_field, i, j, k, grad);
    getVar(i_field, i, j, k-1, var);
  } else {
    getVar(i_field, i, j, k+1, grad);
    getVar(i_field, i, j, k, var);
  }
  for (int i_var = 0; i_var < numVariables(); ++i_var) {
    grad[i_var] -= var[i_var];
    grad[i_var] *= D;
  }
  countFlops(2*numVariables());
  delete [] var;
}

inline bool CartesianPatch::checkRange(size_t i, size_t j, size_t k)
{
  GlobalDebug::ijk(i, j, k);
  if (i >= sizeI() || j >= sizeJ() || k >= sizeK()) {
    return false;
  }
  return true;
}

#ifdef WITH_VTK
template <typename TVariables>
void CartesianPatch::writeToVtk(size_t i_field, QString file_name, const TVariables& variables)
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

  real* raw_var = new real [numVariables()];
  for (int i_var = 0; i_var < variables.numScalars(); ++i_var) {
    vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
    var->SetName(qPrintable(variables.getScalarName(i_var)));
    var->SetNumberOfValues(variableSize());
    grid->GetCellData()->AddArray(var);
    vtkIdType id = 0;
    for (size_t k = 0; k < m_NumK; ++k) {
      for (size_t j = 0; j < m_NumJ; ++j) {
        for (size_t i = 0; i < m_NumI; ++i) {
          getVar(i_field, i, j, k, raw_var);
          var->SetValue(id, variables.getScalar(i_var, raw_var));
          ++id;
        }
      }
    }
  }
  for (int i_var = 0; i_var < variables.numVectors(); ++i_var) {
    vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
    var->SetName(qPrintable(variables.getVectorName(i_var)));
    var->SetNumberOfComponents(3);
    var->SetNumberOfTuples(variableSize());
    grid->GetCellData()->AddArray(var);
    vtkIdType id = 0;
    for (size_t k = 0; k < m_NumK; ++k) {
      for (size_t j = 0; j < m_NumJ; ++j) {
        for (size_t i = 0; i < m_NumI; ++i) {
          getVar(i_field, i, j, k, raw_var);
          vec3_t v = variables.getVector(i_var, raw_var);
          float vf[3];
          vf[0] = v[0]; vf[1] = v[1]; vf[2] = v[2];
          var->SetTuple(id, vf);
          ++id;
        }
      }
    }
  }
  delete [] raw_var;

  vtkSmartPointer<vtkXMLRectilinearGridWriter> vtr = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
  vtr->SetFileName(qPrintable(file_name + ".vtr"));
  vtr->SetInput(grid);
  vtr->Write();
}
#endif


#endif // CARTESIANPATCH_H
