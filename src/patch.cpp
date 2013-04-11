#include "patch.h"

Patch::Patch(size_t num_protectlayers, size_t num_overlaplayers)
{
  m_Data = NULL;
  m_NumFields = 0;
  m_NumVariables = 0;
  m_VariableSize = 0;
  m_FieldSize = 0;
  m_ProtectException = false;
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

bool Patch::readFromFile(istringstream iss_input)
{
  vec3_t xyzoref;        // reference point in parental coords
  vec3_t base_i, base_j; // base vectors of bloc orientation in parental coords.
  iss_input >> xyzoref[0];
  iss_input >> xyzoref[1];
  iss_input >> xyzoref[2];
  iss_input >> base_i[0];
  iss_input >> base_i[1];
  iss_input >> base_i[2];
  iss_input >> base_j[0];
  iss_input >> base_j[1];
  iss_input >> base_j[2];
  iss_input >> m_ioscale;
  setupTransformation(xyzoref,
                      base_i, base_j);
  /// @todo check before returning "true"
  return true;
}

void Patch::setupTransformation(vec3_t xyzoref,
                                vec3_t base_i, vec3_t base_j)
{
  // Build up transformation matrix
  //.. build preliminary transformation this --> inertial
  //   transforms point or freevec given in coord-syst of "this" to
  //   parental coord-syst.
  /// @todo needs testing
  CoordTransform this2inertial;
  this2inertial.setVector(xyzoref);
  this2inertial.setMatrixFromBaseIJ(base_i, base_j);
  //.. build transformation matrix inertial --> this
  m_transformInertial2This.setAll(this2inertial.inverse());

  Transformation t;    /// @todo for compatibility only. Get rid if no longer needed
  t.setVector(xyzoref);
  setTransformation(t.inverse());
}
//void Patch::insertNeighbour(Patch* neighbour_patch) {
//  // Check, if boundary cells are yet extracted. If not, do so.
//  if(!m_receiveCells_OK) {
//    m_receive_cells.clear();
//    extractReceiveCells();
//    compactReceiveCellLists();
//    m_receiveCells_OK = true;
//  }
//  pair<Patch*, CoordTransformVV> dependency_h;
//  dependency_h.first = neighbour_patch;
//  dependency_h.second.setTransFromTo(m_transformInertial2This, neighbour_patch->m_transformInertial2This);
//  m_neighbours.push_back(dependency_h);
//  //  if(m_InterpolateData) {
//  //    InterCoeffWS icws_h;
//  //    m_InterCoeffData_WS.push_back(icws_h);   /// @todo would be nice to have a link_vector<T>
//  //  }
//  //  if(m_InterpolateGrad1N) {
//  //    InterCoeffWS icws_h;
//  //    m_InterCoeffGrad1N_WS.push_back(icws_h);
//  //  }
//  if(m_InterpolateData) {m_InterCoeffData_WS.resize(m_neighbours.size());}   /// @todo would be nice to have a link_vector<T>
//  if(m_InterpolateGrad1N) {m_InterCoeffGrad1N_WS.resize(m_neighbours.size());}
//  size_t i_neighbour = m_neighbours.size() - 1;  /// @todo any better idea to get direct pos index?
//  computeDependencies(i_neighbour);
//}

// changed 2012_08_10: ommit m_neighbours and store in
void Patch::insertNeighbour(Patch* neighbour_patch) {
  // Check, if boundary cells are yet extracted. If not, do so.
  if(!m_receiveCells_OK) {
    m_receive_cells.clear();
    extractReceiveCells();
    compactReceiveCellLists();
    m_receiveCells_OK = true;
  }
  // Insert donor patch and relative coord transformation into m_neighbours and get neighbourship index
  pair<Patch*, CoordTransformVV> dependency_h;
  dependency_h.first = neighbour_patch;
  CoordTransformVV ctvv_relative;
  ctvv_relative.setTransFromTo(m_transformInertial2This, neighbour_patch->m_transformInertial2This);
  dependency_h.second = ctvv_relative;
  m_neighbours.push_back(dependency_h);
  size_t i_neighbour = m_neighbours.size() - 1;  // checked: OK, since size() >= 1

  // Depending on calculation type, insert donor and relative coord transformation into InterCoeffWS lists
  if (m_InterpolateData) {
    m_InterCoeffData_WS.resize(i_neighbour + 1);
    //CHANGE_DONOR! m_InterCoeffData_WS[i_neighbour].setDonorPatch(neighbour_patch);
    m_InterCoeffData_WS[i_neighbour].setCoordTransform(ctvv_relative.extractReverse());
  }   /// @todo would be nice to have a link_vector<T>

  if (m_InterpolateGrad1N) {
    m_InterCoeffGrad1N_WS.resize(i_neighbour + 1);
    //CHANGE_DONOR! m_InterCoeffGrad1N_WS[i_neighbour].setDonorPatch(neighbour_patch);
    m_InterCoeffGrad1N_WS[i_neighbour].setCoordTransform(ctvv_relative.extractReverse());
  }

  bool found_dependency = computeDependencies(i_neighbour);

  // discard donor neighbourship relation, if not any was found
  if (!found_dependency) {
    m_neighbours.pop_back();
    if (m_InterpolateGrad1N) {
      m_InterCoeffGrad1N_WS.pop_back();
    }
    if (m_InterpolateData) {
      m_InterCoeffData_WS.pop_back();
    }
  }
  /** @todo Must apply InterCoeffxxxx_WS::adjust2Average AFTER all neighbouring donor patches have been inserted.
    *       This ensures, the correct number of hits (in case of multiple donor patches contributing to single cells)
    *       is known. */  // DONE!!!
}

void Patch::finalizeDependencies()
{
  // Eliminate pretended receiving cells, if these received no hits at all,
  // m_receive_cell_data_hits[i] == 0 && m_receive_cell_grad1N_hits[i] == 0
  // Do the following:
  //  1) Build a scratch list v_shift with index shifts
  //  2) Build a scratch list index_new_from_old, so that an operation of type
  //     m_receive_cells..new[ll] = m_receive_cells..old[index_new_from_old[ll]]
  //     eliminates all zeros
  //  3) Apply on m_receive_cell_data_hits and m_receive_cell_grad1N_hits
  //  4) Apply on indirect receive indicees stored in InterCoeffWS

  //    ATTENTION: Actually indirect_receiveindex is needed only until InterCoeffWS::adjust2Average is done

  //  1) Build a scratch list v_shift with index shifts
  size_t shift = 0;
  vector<size_t> v_shift;
  v_shift.resize(m_receive_cells.size(), 0);
  for(size_t i_rec = 0; i_rec < m_receive_cells.size(); i_rec++) {
    v_shift[i_rec] = shift;
    bool any_hit = false;
    if(m_InterpolateData) {
      if(m_receive_cell_data_hits[i_rec] > 0) {
        any_hit = true;
      }
    }
    if(m_InterpolateGrad1N) {
      if(m_receive_cell_grad1N_hits[i_rec] > 0) {
        any_hit = true;
      }
    }
    if(!any_hit) { // not any donor contribution found
      shift++;
    }
  }
  //  2) Build a scratch list index_new_from_old, so that an operation of type
  //     m_receive_cells..new[ll] = m_receive_cells..old[index_new_from_old[ll]]
  //     eliminates all zeros
  vector<size_t> index_new_from_old;
  index_new_from_old.resize(m_receive_cells.size(), 0);
  for(size_t i_rec = 0; i_rec < m_receive_cells.size(); i_rec++) {
    index_new_from_old[i_rec] = i_rec - v_shift[i_rec];
  }
  //.. Total number of cells with neighbour contribs
  size_t num_cells_w_contrib = m_receive_cells.size() - v_shift[m_receive_cells.size() - 1];

  //  3) Apply on m_receive_cell_data_hits and m_receive_cell_grad1N_hits
  for(size_t i_rec = 0; i_rec < m_receive_cells.size(); i_rec++) {
    m_receive_cells[index_new_from_old[i_rec]] = m_receive_cells[i_rec];
    if(m_InterpolateData) {
      m_receive_cell_data_hits[index_new_from_old[i_rec]] = m_receive_cell_data_hits[i_rec];
    }
    if(m_InterpolateGrad1N) {
      m_receive_cell_grad1N_hits[index_new_from_old[i_rec]] = m_receive_cell_grad1N_hits[i_rec];
    }
  }
  //  3) Apply on m_receive_cell_data_hits and m_receive_cell_grad1N_hits
  m_receive_cells.resize(num_cells_w_contrib);
  if(m_InterpolateData) {
      m_receive_cell_data_hits.resize(num_cells_w_contrib);
  }
  if(m_InterpolateGrad1N) {
      m_receive_cell_grad1N_hits.resize(num_cells_w_contrib);
  }

//  Hier gehts weiter


  // Reduce contribution weights for receiving cells, being influenced by more than one donor patch.
  for (size_t i_donor = 0; i_donor < m_neighbours.size(); i_donor++) {
    if(m_InterpolateData) {
      m_InterCoeffData_WS[i_donor].reindexIndReceive(index_new_from_old);
      m_InterCoeffData_WS[i_donor].adjust2Average(m_receive_cell_data_hits);
    }
    if(m_InterpolateGrad1N) {
      m_InterCoeffGrad1N_WS[i_donor].reindexIndReceive(index_new_from_old);
      m_InterCoeffGrad1N_WS[i_donor].adjust2Average(m_receive_cell_grad1N_hits);
    }
  }
  // Transfer data to padded data sets, if required.
  if(m_TransferPadded) {    /// @todo I dislike the if-trees. Too many options. Get experience and discard options.
    if(m_InterpolateData) {
      m_InterCoeffData.resize(m_InterCoeffData_WS.size());
      for (size_t i_donor = 0; i_donor < m_neighbours.size(); i_donor++) {
        m_InterCoeffData[i_donor].build(m_InterCoeffData_WS[i_donor]);
        // postponed delete m_InterCoeffData_WS[i_donor];  // erase previous WS-version to save some memory on CPU.
      }
    }
    if(m_InterpolateGrad1N) {
      m_InterCoeffGrad1N.resize(m_InterCoeffGrad1N_WS.size());
      for (size_t i_donor = 0; i_donor < m_neighbours.size(); i_donor++) {
        m_InterCoeffGrad1N[i_donor].build(m_InterCoeffGrad1N_WS[i_donor]);
        // postponed delete m_InterCoeffGrad1N_WS[i_donor];  // erase previous WS-version to save some memory on CPU.
      }
    }
  }
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

/// @todo new mem-structure: must change access to borrowed pointers
void Patch::accessDonorData_WS(const size_t& field)
{
  vector<real*> donor_vars;
  vector<real*> this_vars;
  // assign variable pointers to work on
  for(size_t i_v=0; i_v<numVariables(); i_v++) {
    this_vars.push_back(getVariable(field, i_v));
  }
  // set all receiving data variables to 0, as donors will add their contributions onto
  for(size_t ll_rc=0; ll_rc < m_receive_cells.size(); ll_rc++) {
    /// @todo Eliminate cells from m_receive_cells that did not find any donor neighbour. Eliminates next line.
    if(m_receive_cell_data_hits[ll_rc] != 0) { // only if at least one donor was found. Else do not alter variables.
      size_t l_rc = m_receive_cells[ll_rc];  /// @todo indirect operation!!!
      for(size_t i_v=0; i_v<numVariables(); i_v++) {
        this_vars[i_v][l_rc] = 0.;
      }
    }
  }
  // loop through neighbouring donor patches
  for(size_t i_pd=0; i_pd<m_neighbours.size(); i_pd++) {
    Patch* donor = m_neighbours[i_pd].first;
    //.. assign foreign variable pointers (same sequence as in "this_vars"
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      donor_vars.push_back(donor->getVariable(field, i_v));
    }
    //.. execute transfer contribution
    InterCoeffWS icws = m_InterCoeffData_WS[i_pd];
    icws.transFromTo(donor_vars, this_vars);    // NOTE: NO turning of vectorial variables
  }
}


void Patch::accessTurnDonorData_WS(const size_t& field,
                                   const size_t& i_vx, const size_t& i_vy, const size_t& i_vz)
{
  /// @todo must check performance issues of these little vector allocations, assume OK
  vector<real*> donor_vars;
  vector<real*> this_vars;
  // sort variable sequences to get vectorial quantity at position 0,1,2 in sequence of donor_vars/this_vars
  vector<bool> blocked(numVariables(), false);
  vector<size_t> ind_iv;
  ind_iv.push_back(i_vx);
  ind_iv.push_back(i_vy);
  ind_iv.push_back(i_vz);
  blocked[i_vx] = true;
  blocked[i_vy] = true;
  blocked[i_vz] = true;
  for (size_t i_v = 0; i_v < numVariables(); i_v++) {
    if (!blocked[i_v]) {
      ind_iv.push_back(i_v);
    }
  }
  // assign variable pointers to work on
  for(size_t ii_v=0; ii_v<numVariables(); ii_v++) {
    size_t i_v = ind_iv[ii_v];
    this_vars.push_back(getVariable(field, i_v));
  }
  // set all receiving data variables to 0, as donors will add their contributions onto
  for(size_t ll_rc=0; ll_rc < m_receive_cells.size(); ll_rc++) {
    size_t l_rc = m_receive_cells[ll_rc];
    for(size_t i_v=0; i_v<numVariables(); i_v++) {
      this_vars[i_v][l_rc] = 0.;
    }
  }
  // loop through neighbouring donor patches
  for(size_t i_pd=0; i_pd<m_neighbours.size(); i_pd++) {
    Patch* donor = m_neighbours[i_pd].first;
    //.. assign foreign variable pointers (same sequence as in "this_vars"
    for(size_t ii_v=0; ii_v<numVariables(); ii_v++) {
      size_t i_v = ind_iv[ii_v];
      donor_vars.push_back(donor->getVariable(field, i_v));
    }
    //.. execute transfer contribution
    InterCoeffWS icws = m_InterCoeffData_WS[i_pd];
    icws.transTurnFromTo(donor_vars, this_vars);
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
