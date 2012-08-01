#include "cartesianpatch.h"

#ifdef WITH_VTK
#include <vtkSmartPointer.h>
#include <vtkRectilinearGrid.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <QVector>
#endif

CartesianPatch::CartesianPatch()
{
  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
}

void CartesianPatch::computeDeltas()
{
  double Lx = sqrt(sqr(m_Uxo) + sqr(m_Uyo) + sqr(m_Uzo));
  double Ly = sqrt(sqr(m_Vxo) + sqr(m_Vyo) + sqr(m_Vzo));
  double Lz = sqrt(sqr(m_Wxo) + sqr(m_Wyo) + sqr(m_Wzo));
  countFlops(15);
  countSqrts(3);

  m_Dx = Lx/m_NumI;
  m_Dy = Ly/m_NumJ;
  m_Dz = Lz/m_NumK;
  countFlops(3);

  m_InvDx = 1.0/m_Dx;
  m_InvDy = 1.0/m_Dx;
  m_InvDz = 1.0/m_Dx;
  countFlops(3);

  m_Dixo = m_Uxo/m_NumI;
  m_Diyo = m_Uyo/m_NumI;
  m_Dizo = m_Uzo/m_NumI;
  countFlops(3);

  m_Djxo = m_Vxo/m_NumJ;
  m_Djyo = m_Vyo/m_NumJ;
  m_Djzo = m_Vzo/m_NumJ;
  countFlops(3);

  m_Djxo = m_Wxo/m_NumK;
  m_Djyo = m_Wyo/m_NumK;
  m_Djzo = m_Wzo/m_NumK;
  countFlops(3);

  m_Xco = m_Xo + 0.5*(m_Dixo + m_Djxo + m_Dkxo);
  m_Yco = m_Yo + 0.5*(m_Diyo + m_Djyo + m_Dkyo);
  m_Zco = m_Zo + 0.5*(m_Dizo + m_Djzo + m_Dkzo);
  countFlops(12);
}

void CartesianPatch::setupAligned(real x1, real y1, real z1, real x2, real y2, real z2)
{
  m_Xo = x1;
  m_Yo = y1;
  m_Zo = z1;

  m_Uxo = x2 - x1;
  m_Uyo = 0;
  m_Uzo = 0;
  countFlops(1);

  m_Vxo = 0;
  m_Vyo = y2 - y1;
  m_Vzo = 0;
  countFlops(1);

  m_Wxo = 0;
  m_Wyo = 0;
  m_Wzo = z2 - z1;
  countFlops(1);

  computeDeltas();

  Transformation t;
  t.setVector(vec3_t(x1, y1, z1));
  setTransformation(t.inverse());
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
            var->SetValue(id, f(i_field, i_var, i, j, k));
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
