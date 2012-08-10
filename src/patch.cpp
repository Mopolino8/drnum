#include "patch.h"

Patch::Patch(size_t num_protectlayers, size_t num_overlaplayers)
{
  m_Data = NULL;
  m_NumFields = 0;
  m_NumVariables = 0;
  m_VariableSize = 0;
  m_FieldSize = 0;
  m_NumProtectLayers = num_protectlayers;
  m_NumOverlapLayers = num_overlaplayers;
  m_InterpolateData = false;
  m_InterpolateGrad1N = false;
  m_receiveCells_OK = false;
  m_bbox_OK = false;
}

Patch::~Patch()
{
  deleteData();
}

void Patch::insertNeighbour(Patch* neighbour_patch) {
  // Check, if boundary cells are yet extracted. If not, do so.
  if(!m_receiveCells_OK) {
    m_receive_cells.clear();
    extractReceiveCells();
    compactReceiveCellLists();
    m_receiveCells_OK = true;
  }
  pair<Patch*, CoordTransformVV> dependency_h;
  dependency_h.first = neighbour_patch;
  dependency_h.second.setTransFromTo(m_transformInertial2This, neighbour_patch->m_transformInertial2This);
  m_neighbours.push_back(dependency_h);
//  if(m_InterpolateData) {
//    InterCoeffWS icws_h;
//    m_InterCoeffData_WS.push_back(icws_h);   /// @todo would be nice to have a link_vector<T>
//  }
//  if(m_InterpolateGrad1N) {
//    InterCoeffWS icws_h;
//    m_InterCoeffGrad1N_WS.push_back(icws_h);
//  }
  if(m_InterpolateData) {m_InterCoeffData_WS.resize(m_neighbours.size());}   /// @todo would be nice to have a link_vector<T>
  if(m_InterpolateGrad1N) {m_InterCoeffGrad1N_WS.resize(m_neighbours.size());}
  size_t i_neighbour = m_neighbours.size() - 1;  /// @todo any better idea to get direct pos index?
  computeDependencies(i_neighbour);
}

void Patch::compactReceiveCellLists()
{
  // Clean up m_receive_cells and hit counter lists
  // Sort cell indicees
  // Note: indirect indexing, but piecewise constant strides depending on patch type
  //       example structured: strides with (m_NumI*m_NumJ or m_NumJ*m_NumK or m_NumK*m_NumI)
  sort(m_receive_cells.begin(), m_receive_cells.end());
  // Remove duplicates
  vector<size_t>::iterator it;
  it = unique(m_receive_cells.begin(), m_receive_cells.end());
  m_receive_cells.resize(it - m_receive_cells.begin());
  // Create hit counter lists with same size
  if(m_InterpolateData) {
    m_receive_cell_data_hits.resize(m_receive_cells.size(), 0);
  }
  if(m_InterpolateGrad1N) {
    m_receive_cell_grad1N_hits.resize(m_receive_cells.size(), 0);
  }
}

void Patch::allocateData()
{
  m_FieldSize = m_NumVariables * m_VariableSize;
  m_Data = new real [m_NumFields*m_FieldSize];
}

void Patch::deleteData()
{
  if(m_Data) {   /// @todo condition needed ??
    delete [] m_Data;
  }
}

void Patch::resize(size_t variable_size)
{
  m_VariableSize = variable_size;
  allocateData();
}

void Patch::computeVariableDifference(size_t i_field1, size_t i_var1, size_t i_field2, size_t i_var2, real &max_norm, real &l2_norm)
{
  RESTRICT real *var1 = getVariable(i_field1, i_var1);
  RESTRICT real *var2 = getVariable(i_field2, i_var2);
  max_norm = 0.0;
  l2_norm  = 0.0;
  for (size_t i = 0; i < m_VariableSize; ++i) {
    real diff = sqr(var2[i] - var1[i]);
    max_norm = max(max_norm, diff);
    l2_norm += diff;
  }
  max_norm = sqrt(max_norm);
  l2_norm  = sqrt(l2_norm);
}
