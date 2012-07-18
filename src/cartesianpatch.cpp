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

void CartesianPatch::allocateData(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  Patch::allocateData(m_NumI*m_NumJ*m_NumK);
}

void CartesianPatch::computeDeltas()
{
  double Lx = sqrt(sqr(m_UX) + sqr(m_UY) + sqr(m_UZ));
  double Ly = sqrt(sqr(m_VX) + sqr(m_VY) + sqr(m_VZ));
  double Lz = sqrt(sqr(m_WX) + sqr(m_WY) + sqr(m_WZ));

  m_DX = Lx/m_NumI;
  m_DY = Ly/m_NumJ;
  m_DZ = Lz/m_NumK;

  m_InvDX = 1.0/m_DX;
  m_InvDY = 1.0/m_DX;
  m_InvDZ = 1.0/m_DX;
}

void CartesianPatch::setupAligned(real x1, real y1, real z1, real x2, real y2, real z2)
{
  m_X0 = x1;
  m_Y0 = y1;
  m_Z0 = z1;

  m_UX = x2 - x1;
  m_UY = 0;
  m_UZ = 0;

  m_VX = 0;
  m_VY = y2 - y1;
  m_VZ = 0;

  m_WX = 0;
  m_WY = 0;
  m_WZ = z2 - z1;

  computeDeltas();
}

void CartesianPatch::resize(size_t num_i, size_t num_j, size_t num_k)
{
  /// @todo implement a field mapping mechanism if the grid gets resized
  deleteData();
  allocateData(num_i, num_j, num_k);
  computeDeltas();
}

#ifdef WITH_VTK

void CartesianPatch::writeToVtk(QString file_name)
{
  vtkSmartPointer<vtkRectilinearGrid> grid = vtkSmartPointer<vtkRectilinearGrid>::New();
  QVector<float> x(m_NumI + 1);
  QVector<float> y(m_NumI + 1);
  QVector<float> z(m_NumI + 1);

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
    vtkSmartPointer<vtkFloatArray> var = vtkSmartPointer<vtkFloatArray>::New();
    QString field_name;
    field_name.setNum(i_field + 1);
    if (field_name.size() < 2) {
      field_name = "0" + field_name;
    }
    field_name = "q_" + field_name;
    var->SetName(qPrintable(field_name));
    var->SetNumberOfValues(fieldSize());
    grid->GetCellData()->AddArray(var);
    vtkIdType i_var = 0;
    for (size_t k = 0; k < m_NumK; ++k) {
      for (size_t j = 0; j < m_NumJ; ++j) {
        for (size_t i = 0; i < m_NumI; ++i) {
          var->SetValue(i_var, f(getField(i_field), i, j, k));
          ++i_var;
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
