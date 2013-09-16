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
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "blockobject.h"


BlockObject::BlockObject(PatchGrid* patch_grid,
                         size_t num_layers, size_t num_inner,
                         size_t grey_resolution)
{
  m_PatchGrid = patch_grid;
  m_NumBlackLayers = num_layers;
  m_NumBlackInner = num_inner;
  m_GreyResolution = grey_resolution;
  if (grey_resolution > 20) {
    cout << " Warning, large grey resolution value "
         << grey_resolution
         << " may decay set up performance"
         << endl;
  }
  /// @todo need handling for field index. Allways 0 ?
  m_FieldIndex = 0;
}


//void BlockObject::copyToHost()
//{
//  for (size_t i = 0; i < m_AffectedPatchIDs.size(); ++i) {
//    m_PatchGrid->getPatch(m_AffectedPatchIDs[i])->copyFieldToHost(m_FieldIndex);
//  }
//}

//void BlockObject::copyToDevice()
//{
//  for (size_t i = 0; i < m_AffectedPatchIDs.size(); ++i) {
//    m_PatchGrid->getPatch(m_AffectedPatchIDs[i])->copyFieldToDevice(m_FieldIndex);
//  }
//}


//void BlockObject::operator()()
//{
//  //on GPU copyToHost();

//  // ...

//  //on GPU copyToDevice();
//}


void BlockObject::update(ObjectDefinition *object)
{

  m_ObjectDefinition = object;

  /**
    * @todo Clear method to give number of cells in patch is missing.
    * variableSize() ???
    */

  m_AffectedPatchIDs.clear();  // just to be sure

  m_NumBlackGroups = m_NumBlackLayers + m_NumBlackInner; // m_NumBlackInner may rise, if needed later on
  m_BlackCells.resize(m_NumBlackGroups);

  size_t max_marker     = 10000;  // assume larger than any reasonable group index
  size_t pre_light_grey =  3000;  // assume larger than any reasonable group index
  size_t pre_dark_grey  =  6000;  // assume larger than any reasonable group index


  // loop for patches to find cells located in object definition space
  for (size_t i_p = 0; i_p < m_PatchGrid->getNumPatches(); i_p++) {

    Patch* patch = m_PatchGrid->getPatch(i_p);

    // build a scratch size_t list to mark cells affected
    // cell_marker[any_cell] =
    //   0: marks cell as white (fully outside solid object)
    //   1: marks cell as grey (partially inside and outside solid object)
    //   2 .. m_NumBlackLayers+1: (i-2)th surface layer (behind grey cells)
    //   m_NumBlackLayers+2 .. m_NumBlackInner+1 : i-th inner group

    bool any_hit_in_patch = false;
    bool fully_black = true;

    size_t num_cells = patch->variableSize();
    vector<size_t> cell_marker;
    vector<size_t> grey_count;
    vector<size_t> max_grey_count;
    vector<vec3_t> n_vec;
    cell_marker.resize(num_cells, 0);  // 0 codes as outside object
    grey_count.resize(num_cells);
    max_grey_count.resize(num_cells);
    n_vec.resize(num_cells);

    // Loop for cells of patch to find the ones blocked by object
    // Note 1: dont care black or grey
    // Note 2: dont care subcell-resolution. Take only the cells, whose center is inside
    /// @todo parallel omp makes program extremely slow ???
    //#pragma omp parallel for
    for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
      real xc, yc, zc;
      patch->xyzoCell(l_c,
                      xc, yc, zc);
      if (m_ObjectDefinition->isInside(xc, yc, zc)) {
        cell_marker[l_c] = max_marker;  // preliminary setting to last group
        any_hit_in_patch = true;
      }
      else {
        fully_black = false;
      }
    }

    // Any hit?
    if (any_hit_in_patch) {
      //.. Insert patch in the affected list and get indirect index of patch
      m_AffectedPatchIDs.push_back(i_p);
      size_t ii_p = m_AffectedPatchIDs.size() - 1;

      // Fully black?
      if (fully_black) {
        m_FullyBlackPatches.push_back(ii_p);  // Note double indirect indexing via m_AffectedPatchIDs
      }

      // Build layers. Only, if not fully black or if explicitly requested.
      if (!fully_black) { // Note: On fully black patches, no layers are needed, but inner region groups are buildt later.

        //.. Find frontmost, max_marked cells. These are the ones with any white neighbours.
        //   * Alter marker of the frontmost cells to "preliminary dark grey"
        //   * Mark their previous white neighbours as "preliminary light grey"
        vector<size_t> ind_cell_neighbours;
        for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
          if (cell_marker[l_c] == max_marker) {
            //.... check, if a white cell is in neighbourship
            // patch->cellOverNodeNeighbours(l_c,
            //                               ind_cell_neighbours);
            patch->cellOverFaceNeighbours(l_c,
                                          ind_cell_neighbours);
            bool any_light_neighbour = false;
            for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
              size_t l_cn = ind_cell_neighbours[ll_cn];
              if (cell_marker[l_cn] == 0) { // unmarked, still considered fluid
                //...... Mark previous white neighbour as "preliminary light grey"
                cell_marker[l_cn] = pre_light_grey;
                any_light_neighbour = true;
              }
              else if (cell_marker[l_cn] == pre_light_grey) { // marked as light grey by annother frontmost cell yet
                any_light_neighbour = true;
              }
            } // ll_cn
            if (any_light_neighbour) {  // l_c is a frontmost cell
              //.... Alter marker of l_c to "preliminary dark grey"
              cell_marker[l_c] = pre_dark_grey;
            }
          }  // (cell_marker[l_c] == max_marker)
        }    // l_c

        // Until here, cells are marked:
        // 0:              Sure fluid.
        // pre_light_grey: Outermost cell layer of "blocked object" . These cells may be or really grey, or white.
        // pre_dark_grey:  Previous frontmost cells of "object", that are now covered by a layer of pre_light_grey cells.
        //                 "pre_dark_grey"-cells may in fact be or really grey or black (1st. layer of black cells)
        // max_marker:     Cells blocked by object, farther apart from the surface.

        //.. Analyse preliminary light grey cells, and decide wether these are
        //   * white: grey_count == 0  =>  reset cell back to fluid (marker = 0)
        //   * grey:  grey_count >  0  =>  set cell to really grey (marker = 1)
        //#pragma omp parallel for
        for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
          if (cell_marker[l_c] == pre_light_grey) {
            analyseGreyCount(patch, l_c,
                             grey_count[l_c], max_grey_count[l_c],
                             n_vec[l_c]);

            /// @todo prevent from storing grey_counts, and grads in total patch.

            if (grey_count[l_c] > 0) { // truly a grey cell
              cell_marker[l_c] = 1;  // grey marker, also 1st layer
            }
            else { // white cell (fluid)
              cell_marker[l_c] = 0;
            }
          }
        }

        // Marker state until here, cells are marked:
        // 0:             Sure fluid.
        // 1:             Grey cells at front.
        // pre_dark_grey: Not altered since last marker state.
        // max_marker:    Not altered since last marker state.

        //.. In rather rare cases (concave surface edges, etc...) it is possible, that allthough (light) grey, a
        //   cell may not have any white neighbour. Even though their grey_count might be < max_grey_count, reset
        //   these cells to black.
        for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
          if (cell_marker[l_c] == 1) {  // grey
            //.... check, if a white cell is in neighbourship
            // patch->cellOverNodeNeighbours(l_c,
            //                               ind_cell_neighbours);
            patch->cellOverFaceNeighbours(l_c,
                                          ind_cell_neighbours);
            bool any_white_neighbour = false;
            for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
              size_t l_cn = ind_cell_neighbours[ll_cn];
              if (cell_marker[l_cn] == 0) { // fluid
                any_white_neighbour = true;
              }
            } // ll_cn
            if (!any_white_neighbour) {  // l_c has not any white neighbour: "back to black" (RIP, Amy)
              //.... Alter marker of l_c to "max_marker"
              cell_marker[l_c] = max_marker;
            }
          }
        }

        // Marker state until here, cells are marked:
        // 0:             Sure fluid.
        // 1:             Grey cells at front. Surely with at least one white (fluid) neighbour.
        // pre_dark_grey: Not altered since last marker state.
        // max_marker:    Not altered since last marker state.

        //.. Analyse preliminary dark grey cells wether or not they have any white neighbour
        //   * any white neighbour exists:  =>  the cell is really grey (marker = 1)
        //   * no white neighbour exists:   =>  reset cell marker to max_marker (black)
        for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
          if (cell_marker[l_c] == pre_dark_grey) {
            //.... check, if a white cell is in neighbourship
            //patch->cellOverNodeNeighbours(l_c,
            //                              ind_cell_neighbours);
            patch->cellOverFaceNeighbours(l_c,
                                          ind_cell_neighbours);
            bool any_white_neighbour = false;
            for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
              size_t l_cn = ind_cell_neighbours[ll_cn];
              if (cell_marker[l_cn] == 0) { // white cell (fluid)
                any_white_neighbour = true;
                break;
              }
            }
            if (any_white_neighbour) {
              cell_marker[l_c] = 1;  // cell considered grey. Note: no matter, if grey_count reveals fully black.
              analyseGreyCount(patch, l_c,
                               grey_count[l_c], max_grey_count[l_c],
                               n_vec[l_c]);
            }
            else {
              cell_marker[l_c] = max_marker;
            }
          }
        }

        // Marker state until here, cells are marked:
        // 0:             Sure fluid.
        // 1:             Grey cells at front. Surely with at least one white (fluid) neighbour.
        // max_marker:    All other cells (in core region of "object")
        // !!Note: All front cells declared grey, even if 100% black (grey_count == max_grey_count) , rare.

        //.. Recursively find black layers
        //
        for (size_t i_layer = 0; i_layer < m_NumBlackLayers; i_layer++) {
          size_t layer_marker = i_layer + 2;  // 0 (fluid) and 1 (grey) already used
          //.... loop cells
          for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
            if (cell_marker[l_c] == max_marker) {  // still the unspecified solid marker
              //...... check, if a neighbour cell with lower marker exists
              // patch->cellOverNodeNeighbours(l_c,
              //                               ind_cell_neighbours);
              patch->cellOverFaceNeighbours(l_c,
                                            ind_cell_neighbours);
              for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
                size_t l_cn = ind_cell_neighbours[ll_cn];
                if (cell_marker[l_cn] < layer_marker) {  // neighbour in lower level
                  cell_marker[l_c] = layer_marker;
#ifdef DEBUG
                  if (cell_marker[l_cn] != (layer_marker - 1)) { // not permitted
                    BUG;
                  }
#else
                  break;
#endif
                }  // if (cell_marker[l_cn] < layer_marker)
              }    // for (size_t ll_cn ...
            }      // if (cell_marker[l_c] == max_marker)
          }        // for (size_t l_c ...
        }          // for (size_t i_layer ...
      }            // if (!fully_black || do_black_patches)

      // Marker state until here, cells are marked:
      // 0:                    Sure fluid.
      // 1:                    Grey cells at front.
      // 2 .. m_NumBlackLayers + 1: Cells in layers between grey and inner core of "object"
      // max_marker:    All other cells (in inner region of "object")
      // !!Exception: For fully black patches no layers exist.

      // The rest of cells form a black inner region. Build non recursive checkerboard like
      // groups.
      bool more_groups_needed = true;
      size_t inner_groups = 0;
      vector<size_t> ind_cell_neighbours;
      while (more_groups_needed) {
        more_groups_needed = false;
        size_t current_group_marker = m_NumBlackLayers + 2 + inner_groups;
        //.... loop for cells
        for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
          if (cell_marker[l_c] == max_marker) {  // still the unspecified solid marker
            //...... check if none of your face-neighbours has current_group_marker yet
            patch->cellOverFaceNeighbours(l_c,
                                          ind_cell_neighbours);
            bool block_recurrence = false;
            for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
              size_t l_cn = ind_cell_neighbours[ll_cn];
              if (cell_marker[l_cn] == current_group_marker) {  // neighbour in currect group
                block_recurrence = true;
              }
            }
            if (!block_recurrence) {  // its ok to insert in this group
              cell_marker[l_c] = current_group_marker;
            }
            else {  // a recurrence had to be blocked, need more groups to clear
              more_groups_needed = true;
            }
          }
        }
        inner_groups++;
      }

//#ifdef DEBUG
      // DEBUG ONLY: find nr. cells per types
      size_t max_marker_after = 0;
      for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
        if (cell_marker[l_c] > max_marker_after) {max_marker_after = cell_marker[l_c];}
      }
      vector<size_t> marker_count;
      marker_count.resize(max_marker_after+1, 0);
      for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
        marker_count[cell_marker[l_c]]++;
      }
      cout << endl;
      cout << "Ghost cells in patch " << i_p << endl;
      for (size_t i_m = 0; i_m < marker_count.size(); i_m++) {
        if (marker_count[i_m] > 0) {
          cout << "  marker " << i_m << ": " << marker_count[i_m] << endl;
        }
      }
//#endif

      // Manage memory
      //
      // Check, if number of groups is higher
      if (inner_groups > m_NumBlackInner) {
        m_NumBlackInner = inner_groups;
        m_NumBlackGroups = m_NumBlackLayers + m_NumBlackInner;
        // resize m_BlackCells to accept higher number of groups
        m_BlackCells.resize(m_NumBlackGroups);
      }
      //
      // Expand m_GreyCells and m_BlackCells["all groups"], to allow insertion of the
      // patch i_p (with ind. index ii_p), no matter if patch is fully black (no grey
      // cells) or does not require all groups.
      // Reason: keep patch index coherence between m_GreyCells and m_BlackCells["all groups"].
      m_GreyCells.resize(ii_p + 1);
      for (size_t i_group = 0; i_group < m_NumBlackGroups; i_group++) {
        m_BlackCells[i_group].resize(ii_p + 1);
      }

      // Until here, all cells are marked with gruop markers for grey, black layers and
      // black inner region. Grey counts are stored, if applicable. All groups non recursive.
      // Fill in m_GreyCells and m_BlackCells.
      for (size_t l_c = 0; l_c < patch->variableSize(); l_c++) {
        //.. Decide upon cell colour
        if (cell_marker[l_c] == 1) { // this is a grey cell
          //.. build up greycell_t for this cell
          greycell_t grey_cell;
          grey_cell.l_cell = l_c;
          //.. find influencing cells (fluid neighbours)
          // patch->cellOverNodeNeighbours(l_c,
          //                               ind_cell_neighbours);
          patch->cellOverFaceNeighbours(l_c,
                                        ind_cell_neighbours);
          for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
            size_t l_cn = ind_cell_neighbours[ll_cn];
            if (cell_marker[l_cn] == 0) {
              grey_cell.influencing_cells.push_back(l_cn);
            }
          }
          if (grey_cell.influencing_cells.size() == 0) {  // should not happen at all
            BUG;
          }
          else {
            grey_cell.average_quote = 1./real(grey_cell.influencing_cells.size());
          }
          grey_cell.grey_factor = real(grey_count[l_c]) / real(max_grey_count[l_c]);
          grey_cell.n_vec = n_vec[l_c];
          //.. insert
          m_GreyCells[ii_p].push_back(grey_cell);
        }

        else if (cell_marker[l_c] > 1 && cell_marker[l_c] < m_NumBlackLayers + 2) { // this is a black cell in a layer
          //.. group counting
          size_t layer = cell_marker[l_c] - 2 ;
          //.. build up blackcell_t for this cell
          blackcell_t black_cell;
          black_cell.l_cell = l_c;
          //.. find influencing cells (fluid neighbours)
          vector<size_t> ind_cell_neighbours;
          // patch->cellOverNodeNeighbours(l_c,
          //                               ind_cell_neighbours);
          patch->cellOverFaceNeighbours(l_c,
                                        ind_cell_neighbours);
          for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
            size_t l_cn = ind_cell_neighbours[ll_cn];
            if (cell_marker[l_cn] == cell_marker[l_c] - 1) {
              black_cell.influencing_cells.push_back(l_cn);
            }
          }
          if (black_cell.influencing_cells.size() == 0) {  // should not happen at all
            BUG;
          }
          else {
            black_cell.average_quote = 1./real(black_cell.influencing_cells.size());
          }
          //.. insert
          m_BlackCells[layer][ii_p].push_back(black_cell);
        }

        else if (cell_marker[l_c] >= m_NumBlackLayers + 2 && cell_marker[l_c] < m_NumBlackLayers + 2 + m_NumBlackInner) {
          // this is a black cell in the inner region (not im the layers)
          //.. group counting
          size_t group = cell_marker[l_c] - 2 ;
          //.. build up blackcell_t for this cell
          blackcell_t black_cell;
          black_cell.l_cell = l_c;
          //.. find influencing cells (fluid neighbours)
          vector<size_t> ind_cell_neighbours;
          patch->cellOverFaceNeighbours(l_c,
                                        ind_cell_neighbours);
          for (size_t ll_cn = 0; ll_cn < ind_cell_neighbours.size(); ll_cn++) {
            size_t l_cn = ind_cell_neighbours[ll_cn];
            black_cell.influencing_cells.push_back(l_cn);
          }
          black_cell.average_quote = 1./real(black_cell.influencing_cells.size());
          //.. insert
          m_BlackCells[group][ii_p].push_back(black_cell);
        } // if (cell_marker[l_c] ... else if (cell_marker[l_c]
      }   // for (size_t l_c
    }     // if (any_hit_in_patch)
  }       // for (size_t i_p

  // Check, if data construction is non recurrent (as it should)
#ifdef DEBUG
  checkRecurrence();
#endif
}


void BlockObject::checkRecurrence()
{
  // Create log-lists per patch
  //do later
  //  size_t num_affected = m_AffectedPatchIDs.size();
  //  vector<pair<size_T, size_t> > minmax_grey_influences;
  //  vector<vector<pair<size_T, size_t> > > minmax_blacklayer_influences;
  //  vector<vector<pair<size_T, size_t> > > minmax_blackinner_influences;

  // Loop through patches
  for (size_t ii_p = 0; ii_p < m_AffectedPatchIDs.size(); ii_p++) {
    size_t i_p = m_AffectedPatchIDs[ii_p];
    Patch* patch = m_PatchGrid->getPatch(i_p);
    size_t num_cells = patch->variableSize();
    //.. Assign vectors
    //    grey_cells = m_GreyCells[ii_p];
    //    for (size_t i_bg = 0; i_bg < m_NumBlackGroups; i_bg++) {
    //       black_cells[i_bg] = m_BlackCells[i_bg][ii_p];
    //    }

    //.. Build a scratch marker field for all cells
    vector<size_t> scratch_marker;
    scratch_marker.resize(num_cells, 0);
    //.. Go through grey cells and mark grey
    for (size_t i_grey = 0; i_grey < m_GreyCells[ii_p].size(); i_grey++) {
      scratch_marker[m_GreyCells[ii_p][i_grey].l_cell] = 1;
    }
    //.. Go through black cells and mark according to group index
    for (size_t i_bg = 0; i_bg < m_NumBlackGroups; i_bg++) {
      size_t black_marker = i_bg + 2;
      for (size_t i_black = 0; i_black < m_BlackCells[i_bg][ii_p].size(); i_black++) {
        scratch_marker[m_BlackCells[i_bg][ii_p][i_black].l_cell] = black_marker;
      }
    }

    //.. Go through grey cells and analyse data dependencies
    size_t my_marker = 1;
    //.... loop for grey cells (one grey layer)
    for (size_t ii_c = 0; ii_c < m_GreyCells[ii_p].size(); ii_c++) {
      //...... loop for influencing cells of i_grey
      size_t influence_count = 0;
      for (size_t ii_n = 0; ii_n < m_GreyCells[ii_p][ii_c].influencing_cells.size(); ii_n++) {
        size_t cell_n = m_GreyCells[ii_p][ii_c].influencing_cells[ii_n];
        size_t marker_n = scratch_marker[cell_n];
        if (marker_n != my_marker - 1) {
          BUG;
        }
        influence_count++;
      }
      if (influence_count == 0) {
        BUG;
      }
    }

    //.. Go through black cells in all black layers and analyse data dependencies
    for (size_t i_bg = 0; i_bg < m_NumBlackLayers; i_bg++) {
      size_t my_marker = i_bg + 2;
      //.... loop for black cells
      for (size_t ii_c = 0; ii_c < m_BlackCells[i_bg][ii_p].size(); ii_c++) {
        //...... loop for influencing cells of i_c ( = m_BlackCells[i_bg][ii_p][ii_c].l_cell )
        size_t influence_count = 0;
        for (size_t ii_n = 0; ii_n < m_BlackCells[i_bg][ii_p][ii_c].influencing_cells.size(); ii_n++) {
          size_t cell_n = m_BlackCells[i_bg][ii_p][ii_c].influencing_cells[ii_n];
          size_t marker_n = scratch_marker[cell_n];
          if (marker_n != my_marker - 1) {
            BUG;
          }
          influence_count++;
        }
        if (influence_count == 0) {
          BUG;
        }
      }
    }

    //.. Go through black cells in inner black groups and analyse data dependencies
    //   Note: they may be dependent on cells of different other groups
    for (size_t i_bg = m_NumBlackLayers; i_bg < m_NumBlackGroups; i_bg++) {
      size_t my_marker = i_bg + 2;
      //.... loop for black cells in inner group
      for (size_t ii_c = 0; ii_c < m_BlackCells[i_bg][ii_p].size(); ii_c++) {
        //...... loop for influencing cells of i_grey
        size_t influence_count = 0;
        size_t max_marker = 0;
        size_t min_marker = m_NumBlackGroups;
        for (size_t ii_n = 0; ii_n < m_BlackCells[i_bg][ii_p][ii_c].influencing_cells.size(); ii_n++) {
          size_t cell_n = m_BlackCells[i_bg][ii_p][ii_c].influencing_cells[ii_n];
          size_t marker_n = scratch_marker[cell_n];
          if (marker_n == my_marker) {
            BUG;
          }
          max_marker = max(max_marker, marker_n);
          min_marker = min(min_marker, marker_n);
          influence_count++;
        }
        if (influence_count == 0) {
          BUG;
        }
        if (max_marker >= m_NumBlackGroups + 2) {
          BUG;
        }
        if (min_marker < m_NumBlackLayers + 1) {
          BUG;
        }
      }  // for (size_t ii_c
    }    // for (size_t i_bg
  }      // for (size_t ii_p
}


void BlockObject::analyseGreyCount(Patch* patch, size_t l_cell,
                                   size_t& grey_count, size_t& max_grey_count,
                                   vec3_t& no_vec)
{
  /** @todo ATTENTION: Need consequent coord-syst and naming convention.
    * Currently runs in xyzo (syst. of origin) */

  // Center of cell
  real xo_c, yo_c, zo_c;
  patch->xyzoCell(l_cell,
                  xo_c, yo_c, zo_c);

  // Raster subcell resolutions in patches
  vector<vec3_t> xyzo_subcells;
  vec3_t ref_dxyzo;
  patch->xyzoSubCellRaster (l_cell, m_GreyResolution,
                            xyzo_subcells,
                            ref_dxyzo);

  max_grey_count = xyzo_subcells.size();
  grey_count = 0;

  no_vec[0] = 0.;
  no_vec[1] = 0.;
  no_vec[2] = 0.;

  for (size_t i_sub = 0; i_sub < max_grey_count; i_sub++) {
    real xo_sc = xyzo_subcells[i_sub][0];
    real yo_sc = xyzo_subcells[i_sub][1];
    real zo_sc = xyzo_subcells[i_sub][2];

    if (m_ObjectDefinition->isInside(xo_sc, yo_sc, zo_sc)) {
      grey_count++;

      real delta_xo = xo_sc - xo_c;
      real delta_yo = yo_sc - yo_c;
      real delta_zo = zo_sc - zo_c;
      no_vec[0] -= delta_xo;  // nxo -= delta_xo * 1.; (value "1" for a hit)
      no_vec[1] -= delta_yo;  // nyo -= delta_yo * 1.;
      no_vec[2] -= delta_zo;  // nzo -= delta_zo * 1.;
    }
  }

  // Scale coord directions but do NOT normalize n vector.
  /// @todo Debug this once more carefully before using gradients
  /// @todo Recover thoery for this directional norming
  no_vec[0] *= (1. / (ref_dxyzo[0] * ref_dxyzo[0]));
  no_vec[1] *= (1. / (ref_dxyzo[1] * ref_dxyzo[1]));
  no_vec[2] *= (1. / (ref_dxyzo[2] * ref_dxyzo[2]));

  // Normalize carefully (by epsing it).
  // Reason: no_vec might be singular (0,0,0), if hits in above discrete integral are symmetric
  //         around cell center in all 3 coordinate directions (example: fully black or fully
  //         white cell).
  real eps = 1.e-5 * (ref_dxyzo[0] + ref_dxyzo[1] + ref_dxyzo[2]);  // see above: no_vec ~ counts/length
  real q_len = no_vec[0]*no_vec[0] + no_vec[1]*no_vec[1] + no_vec[2]*no_vec[2] + eps*eps;
  real len = sqrt(q_len); // allways > 0.
  no_vec[0] /= len;
  no_vec[1] /= len;
  no_vec[2] /= len;
}


void BlockObject::setLayerIndexToVar(const size_t &i_field,
                                     const size_t &i_variable)
{
  real Eps = 1.e-5;

  // real Eps = 7.;
  // loop for all patches
  for (size_t i_p = 0; i_p < m_PatchGrid->getNumPatches(); i_p++) {
    Patch* patch = m_PatchGrid->getPatch(i_p);
    real* var = patch->getVariable(i_field, i_variable);
    // set all var = Eps (means "outside")
    for (size_t l_cell = 0; l_cell < patch->variableSize(); l_cell++) {
      var[l_cell] = Eps;
    }
  }

  // grey cells
  //.. loop for patches affected
  for (size_t ii_p = 0; ii_p < m_AffectedPatchIDs.size(); ii_p++) {
    size_t i_p = m_AffectedPatchIDs[ii_p];
    Patch* patch = m_PatchGrid->getPatch(i_p);
    real* var = patch->getVariable(i_field, i_variable);

    //.... set to grey value
    vector<greycell_t>& patch_cells_grey = m_GreyCells[ii_p] ;
    for (size_t ll_cell = 0; ll_cell < patch_cells_grey.size(); ll_cell++) {
      size_t l_cell = patch_cells_grey[ll_cell].l_cell;
      var[l_cell] = patch_cells_grey[ll_cell].grey_factor;
    }
  }

  // black cells
  //.. loop for layers and inner groups (together: groups)
  for (size_t i_group = 0; i_group < m_NumBlackGroups; i_group++) {
    //..., loop for patches affected
    for (size_t ii_p = 0; ii_p < m_AffectedPatchIDs.size(); ii_p++) {
      size_t i_p = m_AffectedPatchIDs[ii_p];
      Patch* patch = m_PatchGrid->getPatch(i_p);
      real* var = patch->getVariable(i_field, i_variable);

      //.... set to group index
      vector<blackcell_t>& patch_cells_black = m_BlackCells[i_group][ii_p] ;
      for (size_t ll_cell = 0; ll_cell < patch_cells_black.size(); ll_cell++) {
        size_t l_cell = patch_cells_black[ll_cell].l_cell;
        var[l_cell] = real(i_group + 1); // !!note: 1st black layer: 1, 2nd: 2, etc ...
      }
    }
  }
}
