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
#ifndef RUNGEKUTTAPG1_H
#define RUNGEKUTTAPG1_H

#include "timeintegration.h"
#include "patchgrid.h"
#include <list>
#include <vector>


/** Test version for multiple patches:
  *  - do data exchange before iterating
  */

class RungeKuttaPG1 : public TimeIntegration
{

private: // attributes

  list<real> m_Alpha;

  PatchGrid* m_PatchGrid;

  vector<size_t> m_SyncField;


public:

  RungeKuttaPG1();

  RungeKuttaPG1(PatchGrid* patch_grid);

  void addAlpha(real a) { m_Alpha.push_back(a); }

  void setPatchGrid(PatchGrid* a_PatchGrid) { m_PatchGrid = a_PatchGrid; }

  void addSyncField(size_t f) { m_SyncField.push_back(f); }

  virtual void operator()(real dt);

};

#endif // RUNGEKUTTAPG1_H
