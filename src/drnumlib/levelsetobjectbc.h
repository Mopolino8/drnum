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
#ifndef LEVELSETOBJECTBC_H
#define LEVELSETOBJECTBC_H

class LevelSetObjectBC;

#include "genericoperation.h"
#include "levelsetobject.h"
//#include "LSLayerData.h" /// @todo capital letters?
#include "lslayerdataextrapol.h"

/**
  * Base class for levelset based inner boundary conditions.
  * Names of derived classes will carry the name suffix LSOBC.
  *
  * Attention: In present version, data alignment is not restricted to
  *            affected patches, but considers all in PatchGrid.
  */
class LevelSetObjectBC : public GenericOperation
{

  #include "levelsetobjectbc_common.h"

protected: // attributes
//  size_t m_Field;                   /// the variable field to work on
  LevelSetObject* m_LevelSetObject; /// the LevelSet object to work on
//  size_t m_AbuseField;              /// abused field to avoid recursion
  PatchGrid* m_PatchGrid;           /// the PatchGrid to work on

  /** Number of inner layers on boundaries */
  size_t m_NumInnerLayers;

//  /** Storage of levelset layer data for inner cells
//    * Data alignment (all with actual length, no padding)
//    * 1st dim: patch id (direct; 0, 1, 2, ...)
//    * 2nd dim: layer:
//    *       0 :  has at least one face neighbour with negative G-value (inside)
//    *       1 .. m_NumOuterLayers : furter layers farther outside
//    * 3rd dim: cell indices in layer */
//  LSLayerDataExtrapol* m_InnerCellsLayers;

  /** Patch-Layer starting index in m_InnerCellsLayers.
    * Data alignment (num_ptches * num_inner_layers)
    * 1st dim: patch id (direct: 0, 1, 2, ...)
    * 2nd dim: layer index (direct: 0, 1, 2, ...) */
  size_t** m_InnerCLStart;

  /** Patch-Layer starting index in m_InnerCellsLayers.
    * 1D setting: pointers for all inner layers in patch (num_patches)
    * 1st dim: patch id (direct: 0, 1, 2, ...) */
  size_t* m_InnerCLStartAll;

  /** Number of outer layers on boundaries */
  size_t m_NumOuterLayers;

//  /** Storage of levelset layer data for outer cells
//    * Data alignment (all with actual length, no padding)
//    * 1st dim: patch id (direct; 0, 1, 2, ...)
//    * 2nd dim: layer:
//    *       0 :  has at least one face neighbour with negative G-value (inside)
//    *       1 .. m_NumOuterLayers : furter layers farther outside
//    * 3rd dim: cell indices in layer */
//  LSLayerDataExtrapol* m_OuterCellsLayers;

  /** Patch-Layer starting index in m_OuterCellsLayers.
    * Data alignment (num_ptches * num_outer_layers)
    * 1st dim: patch id (direct: 0, 1, 2, ...)
    * 2nd dim: layer index (direct: 0, 1, 2, ...) */
  size_t** m_OuterCLStart;

  /** Patch-Layer starting index in m_OuterCellsLayers.
    * 1D setting: pointers for all outer layers in patch (num_patches)
    * 1st dim: patch id (direct: 0, 1, 2, ...) */
  size_t* m_OuterCLStartAll;


public:

  /** Constructor.
    * @param field the variable field to work on
    * @param levelset_object the LevelSet object to work on
    * @param abuse_field abused variable field, needed to avoid recursion. Must be free.
    */
  LevelSetObjectBC (size_t field,
                    LevelSetObject* levelset_object,
                    size_t abuse_field);

  // void setLevelSetObject (LevelSetObject* levelset_object) {m_LevelSetObject = levelset_object;}

  /** Transfer levelset data (cells, values and grads) as given in m_LevelSetObject to
    * array based data sets m_InnerCellsLayers and m_OuterCellsLayers
    */
  void transferCellLayerData ();

  virtual void operator ()() = 0;

};

#endif // LEVELSETOBJECTBC_H
