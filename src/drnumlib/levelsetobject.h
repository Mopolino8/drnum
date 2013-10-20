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
#ifndef LEVELSETOBJECT_H
#define LEVELSETOBJECT_H

class LevelSetObject;

#include "levelsetdefinition.h"
#include "patchgrid.h"
#include "lslayerdataextrapol.h"

/**
  * Class to hold levelset based object definitions.
  * Present implementation relys on a numerical variable.
  *
  * Data settings are taken over from a corresponding LevelSetDefinition
  *
  * Boundaries of objects are located at G = 0.
  *
  * Outer domain: G > 0. where the value of G corresponds to a distance
  * from the surface as defined by the corresponding LevelSetDefinition
  *
  * Inner domain: G < 0. where the value of G corresponds to a distance
  * from the surface as defined by the corresponding LevelSetDefinition
  */
class LevelSetObject
{
  /** @todo Swith to a leaner data concept, rather than using a numerical variable. */
  //  friend class LevelSetObjectBC;
  //  friend class CompressibleEulerLSOBC;

protected: // attributes

  LevelSetDefinition* m_LevelSetDefinition; /// Geometric levelset definition of object(s)
  PatchGrid* m_PatchGrid;                   /// PatchGrid to work on
  size_t     m_FieldIndex;                  /// Index of the variable field to write data to
  size_t     m_VarIndex;                    /// Index of the variable to write data to
  /** @todo Swith to a leaner data concept, rather than using a numerical variable. */

  vector<size_t> m_AffectedPatchIDs;        /// Indices of patches affected by this LevelSetObject
  vector<size_t> m_FullyBlackPatchIDs;      /// Indices of patches that are totally inside the object
  /** @todo Might be better to exclude the fully black patches from m_AffectedPatchIDs (?). */

  size_t m_NumInnerLayers;                  /// Number of inner cell layers to extract
  size_t m_NumOuterLayers;                  /// Number of outer cell layers to extract
  real   m_MinInnerRelDist;                 /// Minimum relative inner distance for 0th layer cells
  real   m_MinOuterRelDist;                 /// Minimum relative inner distance for 0th layer cells

  /** Cell layers around boundaries towards the inside of the object.
    * All cells will have negative G-value.
    * mem stucture:
    *   1st dim: layer:
    *       0 :  has at least one face neighbour with positive G-value (outside)
    *       1 .. m_NumInnerLayers : furter layers towards inside
    *   2nd dim: cell indices in layer */
  vector<vector<vector<LSLayerDataExtrapol> > > m_InnerCellsLayers;

  /** Cell layers around boundaries outside of the object.
    * All cells will have positive G-value.
    * mem stucture:
    *   1st dim: layer:
    *       0 :  has at least one face neighbour with negative G-value (inside)
    *       1 .. m_NumOuterLayers : furter layers farther outside
    *   2nd dim: cell indices in layer */
  vector<vector<vector<LSLayerDataExtrapol> > > m_OuterCellsLayers;


public: // methods

  /** Constructor.
    * @param levelset the levelset definition of the (combined) object to use.
    * @param patch_grid PatchGrid to work on
    * @param field_index Index of variable field to store G-value
    * @param var_index Index of variable to store G-value
    * @param num_inner_layers number of inner cell layers (when extracting cell layers)
    * @param num_outer_layers number of outer cell layers (when extracting cell layers)
    * @param min_innerreldist minimum relative inner distance for 0th layer inner cells
    * @param min_outerreldist minimum relative outer distance for 0th layer outer cells
    */
  LevelSetObject(LevelSetDefinition* levelset,
                 PatchGrid* patch_grid,
                 size_t field_index,
                 size_t var_index,
                 size_t num_inner_layers,
                 size_t num_outer_layers,
                 real min_innerreldist,
                 real min_outerreldist);
  /** @todo Must control number of layers compared to the number of interpatch overlap layers.
    Otherwise the layer cells would get lost in overlap.

  /** Build up method. Constructs variable G on patch_grid.
    */
  void update();

  /** Extract cell layers around the boundaries.
    */
  void extractBCellLayers();
  /// @todo extractBCellLayers might better go into LevelSetObjectBC or even derived

  PatchGrid* getPatchGrid() {return m_PatchGrid;}

  vector<size_t> getAffectedPatchIDs() {return m_AffectedPatchIDs;}

  vector<size_t>* getAffectedPatchIDsPtr() {return &(m_AffectedPatchIDs);}

  size_t getNumInnerLayers() {return m_NumInnerLayers;}

  size_t getNumOuterLayers() {return m_NumOuterLayers;}

  vector<vector<vector<LSLayerDataExtrapol> > >& getInnerCellsLayers() {return m_InnerCellsLayers;}

  vector<vector<vector<LSLayerDataExtrapol> > >& getOuterCellsLayers() {return m_OuterCellsLayers;}

  vector<vector<vector<LSLayerDataExtrapol> > >* getInnerCellsLayersPtr() {return &m_InnerCellsLayers;}

  vector<vector<vector<LSLayerDataExtrapol> > >* getOuterCellsLayersPtr() {return &m_OuterCellsLayers;}

};

#endif // LEVELSETOBJECT_H
