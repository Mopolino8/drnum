#include "prismaticlayerpatch.h"

//#ifdef WITH_VTK
//#include <vtkSmartPointer.h>
//#include <vtkRectilinearGrid.h>
//#include <vtkXMLRectilinearGridWriter.h>
//#include <vtkFloatArray.h>
//#include <vtkCellData.h>
//#include <QVector>
//#endif


PrismaticLayerPatch::PrismaticLayerPatch(size_t num_seeklayers, size_t num_addprotectlayers)
  : Patch(num_seeklayers, num_addprotectlayers)
{
  m_mypatchtype = 1101;
  m_Eps = 1.e-5; /// @todo need a better eps-handling.
  m_Max2DNodeIndexOK = false;
}

bool PrismaticLayerPatch::readFromFile(istringstream& iss_input)
{
  // read number of layers
  // read reference faces (2D indexing with index set each)
  // read coordinates (full 3D node indexing)
  findMax2DNodeIndex();
  BUG; // not yet implemented
}


bool PrismaticLayerPatch::writeToFile(ifstream &s_mesh)
{
}


void PrismaticLayerPatch::findMax2DNodeIndex()
{
  vector<size_t> &node_list = m_RefFacesTo2DNodes.accessItemList();
  m_Max2DNodeIndex = 0;
  for(size_t ii = 0; ii < node_list.size(); ii++) {
    if(node_list[ii] > m_Max2DNodeIndex) {
      m_Max2DNodeIndex = node_list[ii];
    }
  }
}


void PrismaticLayerPatch::scaleRefParental(real scfactor)
{
  Patch::scaleRefParental(scfactor);
  BUG; // not implemented
}


void PrismaticLayerPatch::setSeekException(const size_t& i_bc,
                                           const size_t& num_seeklayers)
{
  m_SeekExceptions = true;
  // set the individual protection exception for i_bc
  BUG; // not implemented
}


void PrismaticLayerPatch::buildBoundingBox()
{
  // go through all nodes to find limita upon coordinates
  BUG; // not implemented
}


void PrismaticLayerPatch::extractReceiveCells()
{
  BUG; // not implemented
}


bool PrismaticLayerPatch::computeDependencies(const size_t& i_neighbour)
{
  BUG; // not implemented
}


void PrismaticLayerPatch::setupInterpolators()
{
  // nothing to be done here (??)
  BUG; // not implemented
}


bool PrismaticLayerPatch::computeCCDataInterpolCoeffs(real x, real y, real z,
                                                      WeightedSet<real>& w_set)
{
  BUG; // not implemented
}

bool PrismaticLayerPatch::computeCCGrad1NInterpolCoeffs(real x, real y, real z,
                                                        real nx, real ny, real nz,
                                                        WeightedSet<real>& d_dn)
{
  BUG; // not implemented
}

