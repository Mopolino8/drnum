// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "blockobject.h"


BlockObject::BlockObject(PatchGrid* patch_grid)
{
  m_PatchGrid = patch_grid;
  m_FieldIndex = 0;
}


void BlockObject::copyToHost()
{
  for (size_t i = 0; i < m_AffectedPatchIDs.size(); ++i) {
    m_PatchGrid->getPatch(m_AffectedPatchIDs[i])->copyFieldToHost(m_FieldIndex);
  }
}

void BlockObject::copyToDevice()
{
  for (size_t i = 0; i < m_AffectedPatchIDs.size(); ++i) {
    m_PatchGrid->getPatch(m_AffectedPatchIDs[i])->copyFieldToDevice(m_FieldIndex);
  }
}

void BlockObject::processFront(vector<vector<pair<vector<size_t>, real> > > &front, real &p_average, real &T_average)
{
  int N = 0;
  p_average = 0;
  T_average = 0;
}

void BlockObject::operator()()
{
  copyToHost();

  // ...

  copyToDevice();
}


void BlockObject::update()
{
  /**
    * @todo Method may be rewritten to allow arbitrary numbers of front layers.
    * Doing so, allows higher order discretizations (>=3).
    */
  /**
    * @todo Clear method to give number of cells in patch is missing.
    * variableSize() ???
    */

  m_AffectedPatchIDs.clear();  // just to be sure

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
      m_AffectedPatchIDs.push_back(i_p);
      size_t ii_p = m_AffectedPatchIDs.size() - 1;
      m_CellsFront1.resize(ii_p + 1);
      m_CellsFront2.resize(ii_p + 1);
      m_CellsInside.resize(ii_p + 1);

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
            cell_BC_mask.first.push_back(i_c); // note: [0] for index of cell to apply blocking
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
            cell_BC_mask.first.push_back(i_c); // note: [0] for index of cell to apply blocking
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
          cell_BC_mask.first.push_back(i_c); // note: [0] for index of cell to apply blocking
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


void BlockObject::setLayerIndexToVar (real** vvar)
{
  real Eps = 1.e-5;
  // loop for püatches affected
  for (size_t ii_p = 0; ii_p < m_AffectedPatchIDs.size(); ii_p++) {
    size_t i_p = m_AffectedPatchIDs[ii_p];
    Patch* patch = m_PatchGrid->getPatch(i_p);
    real* var = vvar[i_p];

    // set all var = Eps (means "outside")
    for (size_t l_cell = 0; l_cell < patch->variableSize(); l_cell++) {
      var[l_cell] = Eps;
    }

    // cells in front 1 layer:
    {
      vector<pair<vector<size_t>, real> >& patch_cells_front1 = m_CellsFront1[ii_p];
      for (size_t ll_cell = 0; ll_cell < patch_cells_front1.size(); ll_cell++) {
        size_t l_cell = patch_cells_front1[ll_cell].first[0];
        var[l_cell] = 1.;
      }
    }

    // cells in front 2 layer:
    {
      vector<pair<vector<size_t>, real> >& patch_cells_front2 = m_CellsFront2[ii_p];
      for (size_t ll_cell = 0; ll_cell < patch_cells_front2.size(); ll_cell++) {
        size_t l_cell = patch_cells_front2[ll_cell].first[0];
        var[l_cell] = 2.;
      }
    }

    // cells in inside region of block object
    {
      vector<pair<vector<size_t>, real> >& patch_cells_inside = m_CellsInside[ii_p];
      for (size_t ll_cell = 0; ll_cell < patch_cells_inside.size(); ll_cell++) {
        size_t l_cell = patch_cells_inside[ll_cell].first[0];
        var[l_cell] = 3.;
      }
    }
  }
}
