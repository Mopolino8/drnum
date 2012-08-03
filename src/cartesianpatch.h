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

  size_t m_NumI;
  size_t m_NumJ;
  size_t m_NumK;

  real m_Xo;
  real m_Yo;
  real m_Zo;

  real m_Xco;
  real m_Yco;
  real m_Zco;

  real m_Uxo;
  real m_Uyo;
  real m_Uzo;

  real m_Vxo;
  real m_Vyo;
  real m_Vzo;

  real m_Wxo;
  real m_Wyo;
  real m_Wzo;

  real m_Dx;
  real m_Dy;
  real m_Dz;
  real m_InvDx;
  real m_InvDy;
  real m_InvDz;

  real m_Dixo;
  real m_Diyo;
  real m_Dizo;

  real m_Djxo;
  real m_Djyo;
  real m_Djzo;

  real m_Dkxo;
  real m_Dkyo;
  real m_Dkzo;

  real m_LimiterEpsilon;

protected: // methods

  void computeDeltas();

public: // methods

  CartesianPatch();
  void setupAligned(real x1, real y1, real z1, real x2, real y2, real z2);
  void resize(size_t num_i, size_t num_j, size_t num_k);

  size_t sizeI() { return m_NumI; }
  size_t sizeJ() { return m_NumJ; }
  size_t sizeK() { return m_NumK; }

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first Cartesian index
   * @param j second Cartesian index
   * @param k third Cartesian index
   * @return the index in the one dimensional data field
   */
  size_t index(int i, int j, int k) { return i*m_NumJ*m_NumK + j*m_NumK + k; }

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
   * Check if an (i,j,k) triple is inside the patch.
   * Attention only the upper limit will be checked (unsigned data type).
   * @param i first index
   * @param j second index
   * @param k third index
   * @return true if it is a valid (i,j,k) triple, false otherwise
   */
  bool checkRange(size_t i, size_t j, size_t k);

  real xo(size_t i, size_t j, size_t k) { return m_Xco + i*m_Dixo + j*m_Djxo + k*m_Dkxo; }
  real yo(size_t i, size_t j, size_t k) { return m_Yco + i*m_Diyo + j*m_Djyo + k*m_Dkyo; }
  real zo(size_t i, size_t j, size_t k) { return m_Zco + i*m_Dizo + j*m_Djzo + k*m_Dkzo; }

  real dx() { return m_Dx; }
  real dy() { return m_Dy; }
  real dz() { return m_Dz; }
  real dV() { return m_Dx*m_Dy*m_Dz; }
  real idx() { return m_InvDx; }
  real idy() { return m_InvDy; }
  real idz() { return m_InvDz; }

  real dixo() { return m_Dixo; }
  real diyo() { return m_Diyo; }
  real dizo() { return m_Dizo; }
  real djxo() { return m_Djxo; }
  real djyo() { return m_Djyo; }
  real djzo() { return m_Djzo; }
  real dkxo() { return m_Dkxo; }
  real dkyo() { return m_Dkyo; }
  real dkzo() { return m_Dkzo; }

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
  return getVariable(i_field, i_var)[i*m_NumJ*m_NumK + j*m_NumK + k];
}

inline void CartesianPatch::getVar(size_t i_field, size_t i, size_t j, size_t k, real *var)
{
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    var[i_var] = f(i_field, i_var, i, j, k);
  }
}

inline void CartesianPatch::setVar(size_t i_field, size_t i, size_t j, size_t k, real *var)
{
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    f(i_field, i_var, i, j, k) = var[i_var];
  }
}

inline bool CartesianPatch::checkRange(size_t i, size_t j, size_t k)
{
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
