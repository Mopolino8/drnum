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
#ifndef BLOCKOBJECT_H
#define BLOCKOBJECT_H

class BlockObject;

#include "genericoperation.h"
#include "patchgrid.h"
#include "perfectgas.h"

class BlockObject : public GenericOperation
{

protected: // data

  PerfectGas m_Gas;
  PatchGrid* m_PatchGrid;
  size_t     m_FieldIndex;

  // mem stucture for m_Cells... data;
  //  1st dim: counter index of affected patch
  //  2nd dim: counter index of affected cell in patch
  //  3rd dim (pair::first) : [0]: cell index in patch
  //                          [1]..size() : indices of influencing cells
  //          (pair::second): real (1/(size()-1)
  vector<vector<pair<vector<size_t>, real> > > m_CellsFront1;
  vector<vector<pair<vector<size_t>, real> > > m_CellsFront2;
  vector<vector<pair<vector<size_t>, real> > > m_CellsInside;

  // postponed
  //  1st dim: layer (front0, front, 1, ..., inside) [m_NumLayers + 1]
  //  2nd dim: counter index of affected patch
  //  3rd dim: counter index of affected cell in patch
  //  4th dim (pair::first) : [0]: cell index in patch
  //                          [1]..size() : indices of influencing cells
  //          (pair::second): real (1.(size()-2)
  //  vector<vector<vector<pair<vector<size_t>, real> > > > m_CellsAllLayers;


  vector<size_t> m_AffectedPatchIDs;

protected: // methods

  void copyToHost();
  void copyToDevice();
  void processFront(vector<vector<pair<vector<size_t>, real> > >& front, real& p_average, real& T_average);


public: // methods

  BlockObject(PatchGrid* patch_grid);

  void update();

  virtual bool isInside(const real& xo, const real& yo, const real& zo) = 0;

  bool isInside(vec3_t xyzo) { return isInside (xyzo[0], xyzo[1],xyzo[2]); }

  vector<size_t> getAffectedPatchIDs() { return m_AffectedPatchIDs; }

  vector<size_t>* getAffectedPatchIDsPtr() { return &(m_AffectedPatchIDs); }

  void setLayerIndexToVar (real** var);

  virtual void operator()();
};

#endif // BLOCKOBJECT_H
