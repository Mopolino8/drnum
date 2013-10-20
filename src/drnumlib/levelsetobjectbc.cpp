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
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "levelsetobjectbc.h"


LevelSetObjectBC::LevelSetObjectBC (size_t field, LevelSetObject* levelset_object)
{
  m_Field = field;
  m_LevelSetObject = levelset_object;
  m_NumInnerLayers = levelset_object->getNumInnerLayers();
  m_NumOuterLayers = levelset_object->getNumOuterLayers();
  m_PatchGrid      = levelset_object->getPatchGrid();
  //m_NumPatches     = m_PatchGrid->getNumPatches();
}


void LevelSetObjectBC::transferCellLayerData ()
{
  size_t count;
  size_t count_index;
  size_t num_patches = m_PatchGrid->getNumPatches();

  // Inner
  vector<vector<vector<LSLayerDataExtrapol> > >& inner_cls
      = m_LevelSetObject->getInnerCellsLayers();

  //.. Count total number of data entries in LevelSetObject::m_InnerCellsLayers
  count = 0;
  for (size_t i_p = 0; i_p < num_patches; i_p++) {
    for (size_t i_layer = 0; i_layer < m_NumInnerLayers; i_layer++) {
      count += inner_cls[i_p][i_layer].size();
    }
  }

  //.. Allocate memory for m_InnerCellsLayers and m_InnerCLStart
  m_InnerCellsLayers = new LSLayerDataExtrapol[count];
  m_InnerCLStart = new size_t*[num_patches + 1]; // Note: end index for last patch needed
  for (size_t i_p = 0; i_p <= num_patches; i_p++) {
    m_InnerCLStart[i_p] = new size_t[m_NumInnerLayers + 1]; // Note: end index for last layers needed
  }

  //.. Insert in m_InnerCellsLayers and m_InnerCLStart
  count_index = 0;
  m_InnerCLStart[0][0] = count_index;
  for (size_t i_p = 0; i_p < num_patches; i_p++) {
    vector<vector<LSLayerDataExtrapol> >& inner_cls_p = inner_cls[i_p];
    for (size_t i_layer = 0; i_layer < m_NumInnerLayers; i_layer++) {
      vector<LSLayerDataExtrapol>& inner_cls_p_l = inner_cls_p[i_layer];
      //      m_InnerCLStart[i_p][i_layer] = count_index;
      for (size_t ll_c = 0; ll_c < inner_cls_p_l.size(); ll_c++) {
        m_InnerCellsLayers[count_index] = inner_cls_p_l[ll_c];
        count_index++;
      }
      m_InnerCLStart[i_p][i_layer + 1] = count_index;
    }
    m_InnerCLStart[i_p + 1][0] = count_index;
  }
  for (size_t i_layer = 0; i_layer <= m_NumInnerLayers; i_layer++) {
    m_InnerCLStart[num_patches][i_layer] = count_index;  // note same setting for all layers for num_patches
  }

  // Outer
  vector<vector<vector<LSLayerDataExtrapol> > >& outer_cls
      = m_LevelSetObject->getOuterCellsLayers();

  //.. Count total number of data entries in LevelSetObject::m_OuterCellsLayers
  count = 0;
  for (size_t i_p = 0; i_p < num_patches; i_p++) {
    for (size_t i_layer = 0; i_layer < m_NumOuterLayers; i_layer++) {
      count += outer_cls[i_p][i_layer].size();
    }
  }

  //.. Allocate memory for m_OuterCellsLayers and m_OuterCLStart
  m_OuterCellsLayers = new LSLayerDataExtrapol[count];
  m_OuterCLStart = new size_t*[num_patches + 1]; // Note: end index for last patch needed
  for (size_t i_p = 0; i_p <= num_patches; i_p++) {
    m_OuterCLStart[i_p] = new size_t[m_NumOuterLayers + 1]; // Note: end index for last layers needed
  }

  //.. Insert in m_InnerCellsLayers and m_InnerCLStart
  count_index = 0;
  m_OuterCLStart[0][0] = count_index;
  for (size_t i_p = 0; i_p < num_patches; i_p++) {
    vector<vector<LSLayerDataExtrapol> >& outer_cls_p = outer_cls[i_p];
    for (size_t i_layer = 0; i_layer < m_NumOuterLayers; i_layer++) {
      vector<LSLayerDataExtrapol>& outer_cls_p_l = outer_cls_p[i_layer];
      //      m_OuterCLStart[i_p][i_layer] = count_index;
      for (size_t ll_c = 0; ll_c < outer_cls_p_l.size(); ll_c++) {
        m_OuterCellsLayers[count_index] = outer_cls_p_l[ll_c];
        count_index++;
      }
      m_OuterCLStart[i_p][i_layer + 1] = count_index;
    }
    m_OuterCLStart[i_p + 1][0] = count_index;
  }
  for (size_t i_layer = 0; i_layer <= m_NumOuterLayers; i_layer++) {
    m_OuterCLStart[num_patches][i_layer] = count_index;  // note same setting for all layers for num_patches
  }
}
