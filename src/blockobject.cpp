#include "blockobject.h"


BlockObject::BlockObject()
{
  // ???
}


void BlockObject::operator ()()
{
//  // loop for patches affected
//  for (size_t ii_p = 0; ii_p < m_affectedPatchIDs.size(); ii_p++) {

//    size_t field = 0;

//    size_t i_p = m_affectedPatchIDs[ii_p];
//    Patch* patch = m_PatchGrid->getPatch(i_p);
//    real* rho = patch->getVariable(field, 0);
//    real* rhou = patch->getVariable(field, 1);
//    real* rhov = patch->getVariable(field, 2);
//    real* rhow = patch->getVariable(field, 3);
//    real* rhoE = patch->getVariable(field, 4);

//    //.. loop for front cells of patch
//    for (size_t i_c = 0; i_c < m_CellsFront[i_p].size(); i_c++) {
//      size_t l = m_CellsFront[i_p][i_c];
//      // do BC for front cells of block object
//      applyFrontBC(rho, rhou, rhov, rhow, rhoE, patch, l);
//    }
//    for (size_t i_c = 0; i_c < m_CellsInside[i_p].size(); i_c++) {
//      // do BC for non boundary cells inside block object
//      applyInnerBC(rho, rhou, rhov, rhow, rhoE, patch, l);
//    }
//  }
}


void BlockObject::applyFrontBC(const Patch *patch, const size_t &l)
{
//  size_t field = 0;
//  real* rho = patch->getVariable(field, 0);
//  real* rhou = patch->getVariable(field, 1);
//  real* rhov = patch->getVariable(field, 2);
//  real* rhow = patch->getVariable(field, 3);
//  real* rhoE = patch->getVariable(field, 4);
//  m_Gas::conservativeToPrimitive()



}


void BlockObject::applyInnerBC(const Patch *patch, const size_t &l)
{

}


void BlockObject::fillAll()
{
  m_affectedPatchIDs.clear();  // just to be sure

  // build up first dim of m_Cells... data
  m_CellsFront1.resize(m_PatchGrid->getNumPatches());
  m_CellsFront2.resize(m_PatchGrid->getNumPatches());
  m_CellsInside.resize(m_PatchGrid->getNumPatches());

  // loop for patches to find cells located in object definition space
  for (size_t i_p = 0; i_p < m_PatchGrid->getNumPatches(); i_p++) {
    m_CellsFront1[i_p].clear();  // just to be sure
    m_CellsInside[i_p].clear(); // just to be sure
    //Patch* patch = m_PatchGrid->getPatch(i_p);
    /// @todo need a general function in Patch to access x,y,z
    /// @todo restricted to cartesian at present
    CartesianPatch& patch = *(dynamic_cast<CartesianPatch*>(m_PatchGrid->getPatch(i_p)));

    // build a scratch bool list to mark cells affected
    size_t num_cells = patch.sizeL();
    vector<size_t> cell_marker;
    cell_marker.resize(num_cells, 0);  // 0 codes as outside object

    // loop to find cells affected
    for (size_t i_cell = 0; i_cell < patch.sizeI(); i_cell++) {
      for (size_t j_cell = 0; j_cell < patch.sizeJ(); j_cell++) {
        for (size_t k_cell = 0; k_cell < patch.sizeK(); k_cell++) {
          vec3_t xyzo_cell = patch.xyzoCell(i_cell, j_cell, k_cell);
          if (isInside(xyzo_cell)) {
            size_t l = patch.index(i_cell, j_cell, k_cell);
            cell_marker[l] = 3;  // markes as "body" inside region
          }
        }
      }
    }

    // insert cells affected into "front" and "inside" lists

    // Find front_1 cells: These are cells in object, having at least one fluid neighbour.
    for (size_t i_cell = 0; i_cell < patch.sizeI(); i_cell++) {
      for (size_t j_cell = 0; j_cell < patch.sizeJ(); j_cell++) {
        for (size_t k_cell = 0; k_cell < patch.sizeK(); k_cell++) {
          size_t l = patch.index(i_cell, j_cell, k_cell);
          if (cell_marker[l] == 3) {  // a body cell as found above
            size_t l_n;
            bool error;
            bool inside = true;
            // check neighbours
            // i, j, k-1
            l_n = patch.save_index(i_cell, j_cell, k_cell-1, error);
            if (!error) {
              if (cell_marker[l_n] == 0) {  // a fluid cell
                inside = false;
              }
            }
            // i, j, k+1
            l_n = patch.save_index(i_cell, j_cell, k_cell+1, error);
            if (!error) {
              if (cell_marker[l_n] == 0) {  // a fluid cell
                inside = false;
              }
            }
            // i, j-1, k
            l_n = patch.save_index(i_cell, j_cell-1, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 0) {  // a fluid cell
                inside = false;
              }
            }
            // i, j+1, k
            l_n = patch.save_index(i_cell, j_cell+1, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 0) {  // a fluid cell
                inside = false;
              }
            }
            // i-1, j, k
            l_n = patch.save_index(i_cell-1, j_cell, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 0) {  // a fluid cell
                inside = false;
              }
            }
            // i+1, j, k
            l_n = patch.save_index(i_cell+1, j_cell, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 0) {  // a fluid cell
                inside = false;
              }
            }
            // insert in "front" or "neighbour"
            if (inside) {
              cell_marker[l] = 1;  // this qualifies as "front 1" (first layer)
            }
          }   // cell_marker[l] == 3
        }     // k_cell
      }       // j_cell
    }         // i_cell

    // Find front_2 cells: These are cells in object, having at least one neighbour at front 1.
    for (size_t i_cell = 0; i_cell < patch.sizeI(); i_cell++) {
      for (size_t j_cell = 0; j_cell < patch.sizeJ(); j_cell++) {
        for (size_t k_cell = 0; k_cell < patch.sizeK(); k_cell++) {
          size_t l = patch.index(i_cell, j_cell, k_cell);
          if (cell_marker[l] == 3) {  // a body cell as found above
            size_t l_n;
            bool error;
            bool inside = true;
            // check neighbours
            // i, j, k-1
            l_n = patch.save_index(i_cell, j_cell, k_cell-1, error);
            if (!error) {
              if (cell_marker[l_n] == 1) {  // a fluid cell
                inside = false;
              }
            }
            // i, j, k+1
            l_n = patch.save_index(i_cell, j_cell, k_cell+1, error);
            if (!error) {
              if (cell_marker[l_n] == 1) {  // a fluid cell
                inside = false;
              }
            }
            // i, j-1, k
            l_n = patch.save_index(i_cell, j_cell-1, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 1) {  // a fluid cell
                inside = false;
              }
            }
            // i, j+1, k
            l_n = patch.save_index(i_cell, j_cell+1, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 1) {  // a fluid cell
                inside = false;
              }
            }
            // i-1, j, k
            l_n = patch.save_index(i_cell-1, j_cell, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 1) {  // a fluid cell
                inside = false;
              }
            }
            // i+1, j, k
            l_n = patch.save_index(i_cell+1, j_cell, k_cell, error);
            if (!error) {
              if (cell_marker[l_n] == 1) {  // a fluid cell
                inside = false;
              }
            }
            // insert in "front" or "neighbour"
            if (inside) {
              cell_marker[l] = 2;  // this qualifies as "front 2" (second layer)
            }
          }   // cell_marker[l] == 3
        }     // k_cell
      }       // j_cell
    }         // i_cell

    // insert in m_Cells... lists
    for (size_t i_cell = 0; i_cell < patch.sizeI(); i_cell++) {
      for (size_t j_cell = 0; j_cell < patch.sizeJ(); j_cell++) {
        for (size_t k_cell = 0; k_cell < patch.sizeK(); k_cell++) {
          size_t l = patch.index(i_cell, j_cell, k_cell);
          if (cell_marker[l] == 1) {
            m_CellsFront1[i_p].push_back(l);
            size_t ll = m_CellsFront1[i_p].size()-1;
            size_t l_n;
            pair<vector<size_t>, real> contrib;
            l_n = patch.save_index(i_cell, j_cell, k_cell-1, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell, j_cell, k_cell+1, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell, j_cell-1, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell, j_cell+1, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell-1, j_cell, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell+1, j_cell, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
          }
          else if (cell_marker[l] == 2) {
            m_CellsFront2[i_p].push_back(l);
            size_t ll = m_CellsFront2[i_p].size()-1;
            size_t l_n;
            pair<vector<size_t>, real> contrib;
            l_n = patch.save_index(i_cell, j_cell, k_cell-1, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell, j_cell, k_cell+1, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell, j_cell-1, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell, j_cell+1, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell-1, j_cell, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
            l_n = patch.save_index(i_cell+1, j_cell, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }
          }
          else if (cell_marker[l] == 2) {









            if (k_cell != 0) {
              l_n = patch.index(i_cell, j_cell, k_cell-1);
              if (cell_marker[l_n] == 0) {contrib.first.push_back(l_n);}
            }
            if (k_cell != patch.sizeK()-1) {
              l_n = patch.index(i_cell, j_cell, k_cell+1);
              if (cell_marker[l_n] == 0) {contrib.first.push_back(l_n);}
            }


            l_n = patch.save_index(i_cell, j_cell-1, k_cell, error);
            if (!error && cell_marker[l_n] == 0) {
              contrib.first.push_back(l_n);
            }




              if (cell_marker[l_n] == 0) {contrib.first.push_back(l_n);}
            }
            if (k_cell != patch.sizeK()-1) {
              l_n = patch.index(i_cell, j_cell, k_cell+1);
              if (cell_marker[l_n] == 0) {contrib.first.push_back(l_n);}
            }



          }


        }


    // insert patch affected
    if (m_CellsFront[i_p].size() > 0 || m_CellsInside[i_p].size() > 0) {
      m_affectedPatchIDs.push_back(i_p);
    }

  }           // i_p

}
