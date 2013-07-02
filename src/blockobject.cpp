#include "blockobject.h"


BlockObject::BlockObject()
{
  // ???
}


void BlockObject::operator ()()
{
  // field to work on
  size_t field = 0;

  // loop for patches affected
  for (size_t ii_p = 0; ii_p < m_affectedPatchIDs.size(); ii_p++) {

    //.. patch settings
    size_t i_p = m_affectedPatchIDs[ii_p];
    Patch* patch = m_PatchGrid->getPatch(i_p);
    real* rho = patch->getVariable(field, 0);
    real* rhou = patch->getVariable(field, 1);
    real* rhov = patch->getVariable(field, 2);
    real* rhow = patch->getVariable(field, 3);
    real* rhoE = patch->getVariable(field, 4);

    //.. patch related BC masks
    vector<pair<vector<size_t>, real> >& patch_cells_front1 = m_CellsFront1[ii_p];
    vector<pair<vector<size_t>, real> >& patch_cells_front2 = m_CellsFront2[ii_p];
    vector<pair<vector<size_t>, real> >& patch_cells_inside = m_CellsInside[ii_p];

    //.. loop for front_1 cells of patch
    for (size_t ii_c = 0; ii_c < patch_cells_front1.size(); ii_c++) {
      pair<vector<size_t>, real>& cell_front1_data = patch_cells_front1[ii_c];
      size_t i_cell = cell_front1_data.first[0];
      rho[i_cell] = 0.;
      rhou[i_cell] = 0.;
      rhov[i_cell] = 0.;
      rhow[i_cell] = 0.;
      rhoE[i_cell] = 0.;
      //.... Loop for influencing nodes
      for (size_t ii_cn = 1; ii_cn < cell_front1_data.first.size()-1; ii_cn++) {
        size_t i_cn = cell_front1_data.first[ii_cn];
        rho[i_cell]  += rho[i_cn]  * cell_front1_data.second;  // do on primitive vars later
        rhoE[i_cell] += rhoE[i_cn] * cell_front1_data.second;
      }
    }

    //.. loop for front_2 cells of patch
    for (size_t ii_c = 0; ii_c < patch_cells_front2.size(); ii_c++) {
      pair<vector<size_t>, real>& cell_front2_data = patch_cells_front2[ii_c];
      size_t i_cell = cell_front2_data.first[0];
      rho[i_cell] = 0.;
      rhou[i_cell] = 0.;
      rhov[i_cell] = 0.;
      rhow[i_cell] = 0.;
      rhoE[i_cell] = 0.;
      //.... Loop for influencing nodes
      for (size_t ii_cn = 1; ii_cn < cell_front2_data.first.size()-1; ii_cn++) {
        size_t i_cn = cell_front2_data.first[ii_cn];
        rho[i_cell]  += rho[i_cn]  * cell_front2_data.second;  // do on primitive vars later
        rhoE[i_cell] += rhoE[i_cn] * cell_front2_data.second;
      }
    }

    //.. loop for front_2 cells of patch
    //   ATTENTION: recursive
    for (size_t ii_c = 0; ii_c < patch_cells_inside.size(); ii_c++) {
      pair<vector<size_t>, real>& cell_inside_data = patch_cells_inside[ii_c];
      size_t i_cell = cell_inside_data.first[0];
      rho[i_cell] = 0.;
      rhou[i_cell] = 0.;
      rhov[i_cell] = 0.;
      rhow[i_cell] = 0.;
      rhoE[i_cell] = 0.;
      //.... Loop for influencing nodes
      for (size_t ii_cn = 1; ii_cn < cell_inside_data.first.size()-1; ii_cn++) {
        size_t i_cn = cell_inside_data.first[ii_cn];
        rho[i_cell]  += rho[i_cn]  * cell_inside_data.second;  // do on primitive vars later
        rhoE[i_cell] += rhoE[i_cn] * cell_inside_data.second;
      }
    }
  }
}


void BlockObject::fillAll()
{
  /**
    * @todo Method may be rewritten to allow arbitrary numbers of front layers.
    * Doing so, allows higher order discretizations (>=3).
    */

  m_affectedPatchIDs.clear();  // just to be sure

  // build up first dim of m_Cells... data
  m_CellsFront1.resize(m_PatchGrid->getNumPatches());
  m_CellsFront2.resize(m_PatchGrid->getNumPatches());
  m_CellsInside.resize(m_PatchGrid->getNumPatches());

  // loop for patches to find cells located in object definition space
  for (size_t i_p = 0; i_p < m_PatchGrid->getNumPatches(); i_p++) {

    Patch* patch = m_PatchGrid->getPatch(i_p);

    // build a scratch int list to mark cells affected
    //  0: codes as outside object
    //  1: 1st front layer
    //  2: 2nd front layer
    //  3: inside object

    bool any_hit_in_patch = false;
    size_t num_cells = patch->variableSize();
    vector<size_t> cell_marker;
    cell_marker.resize(num_cells, 0);  // 0 codes as outside object

    // Loop for cells of patch to find the ones blocked by object
    // Dont care front1, front2, inside at present
    for (size_t i_c = 0; i_c < patch->variableSize(); i_c++) {
      real xc, yc, zc;
      patch->xyzoCell(i_c,
                      xc, yc, zc);
      if (isInside(xc, yc, zc)) {
        cell_marker[i_c] = 3;  // preliminary setting: "body inside" region
        any_hit_in_patch = true;
      }
    }

    // Any hit?
    if (any_hit_in_patch) {
      // Insert patch in the affected list
      m_affectedPatchIDs.push_back(i_p);
      size_t ii_p = m_affectedPatchIDs.size();
      m_CellsFront1.resize(ii_p);
      m_CellsFront2.resize(ii_p);
      m_CellsInside.resize(ii_p);

      // Resolve first dimension of lists m_Cells...
      vector<pair<vector<size_t>, real> >& patch_cells_front1 = m_CellsFront1[ii_p];
      vector<pair<vector<size_t>, real> >& patch_cells_front2 = m_CellsFront2[ii_p];
      vector<pair<vector<size_t>, real> >& patch_cells_inside = m_CellsInside[ii_p];
      patch_cells_front1.clear();
      patch_cells_front2.clear();
      patch_cells_inside.clear();

      // Find front_1 cells: These are cells in object, having at least one
      // fluid neighbour cell.
      for (size_t i_c = 0; i_c < patch->variableSize(); i_c++) {
        if (cell_marker[i_c] == 3) {
          vector<size_t> ind_cell_neighbours;
          patch->cellNeighbours(i_c,
                                ind_cell_neighbours);
          //.... check neighbours cell_marker
          bool any_fluid = false;
          for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
            size_t l_cell_neighbour = ind_cell_neighbours[ll_cn];
            if (cell_marker[l_cell_neighbour] == 0) {
              any_fluid = true;
              break;
            }
          }
          if (any_fluid) {
            cell_marker[i_c] = 1; // cell on front_1
            //...... resolve
            pair<vector<size_t>, real> cell_BC_mask;
            cell_BC_mask.first[0] = i_c; // note: [0] for index of cell to apply blocking
            //...... loop again to set fluid neighbours
            for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
              size_t l_cell_neighbour = ind_cell_neighbours[ll_cn];
              if (cell_marker[l_cell_neighbour] == 0) {
                cell_BC_mask.first.push_back(l_cell_neighbour);
              }
            }
            cell_BC_mask.second = 1. / float(cell_BC_mask.first.size()-1);
            patch_cells_front1.push_back(cell_BC_mask);
          }
        }
      }

      // Find front_2 cells: These are cells in object, that are not located in
      // front1 layer, however having at least one neighbour cell in front1 laeyer.
      for (size_t i_c = 0; i_c < patch->variableSize(); i_c++) {
        if (cell_marker[i_c] == 3) {
          vector<size_t> ind_cell_neighbours;
          patch->cellNeighbours(i_c,
                                ind_cell_neighbours);
          //.... check neighbours cell_marker
          bool any_lower_layer = false;
          for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
            size_t l_cell_neighbour = ind_cell_neighbours[ll_cn];
            if (cell_marker[l_cell_neighbour] == 1) {
              any_lower_layer = true;
              break;
            }
          }
          if (any_lower_layer) {
            cell_marker[i_c] = 2; // cell on front_2
            //...... resolve
            pair<vector<size_t>, real> cell_BC_mask;
            cell_BC_mask.first[0] = i_c; // note: [0] for index of cell to apply blocking
            //...... loop again to set lower layer neighbours
            for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
              size_t l_cell_neighbour = ind_cell_neighbours[ll_cn];
              if (cell_marker[l_cell_neighbour] == 1) {
                cell_BC_mask.first.push_back(l_cell_neighbour);
              }
            }
            cell_BC_mask.second = 1. / float(cell_BC_mask.first.size()-1);
            patch_cells_front2.push_back(cell_BC_mask);
          }
        }
      }

      // Insert the inner cells of the blockobject in patch_cells_inside
      for (size_t i_c = 0; i_c < patch->variableSize(); i_c++) {
        if (cell_marker[i_c] == 3) {
          vector<size_t> ind_cell_neighbours;
          patch->cellNeighbours(i_c,
                                ind_cell_neighbours);
          //.... resolve
          pair<vector<size_t>, real> cell_BC_mask;
          cell_BC_mask.first[0] = i_c; // note: [0] for index of cell to apply blocking
          //..... loop to set neighbours (all!)
          for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
            size_t l_cell_neighbour = ind_cell_neighbours[ll_cn];
            cell_BC_mask.first.push_back(l_cell_neighbour);
          }
          cell_BC_mask.second = 1. / float(cell_BC_mask.first.size()-1);
          patch_cells_inside.push_back(cell_BC_mask);
        }
      }
    }   // if (any_hit_in_patch)
  }     // for (size_t i_p = 0; i_p < m_PatchGrid->getNumPatches(); i_p++)
}

