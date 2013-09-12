#ifndef BLOCKOBJECT_H
#define BLOCKOBJECT_H

struct greycell_t;
struct blackcell_t;
class BlockObject;

#include "genericoperation.h"
#include "objectdefinition.h"
#include "patchgrid.h"
#include "perfectgas.h"

struct greycell_t
{
  size_t l_cell;                    // ghost cell index in patch indexing
  vector<size_t> influencing_cells; // cells influencing ghost_cell (on same patch only)
  real average_quote;               // 1./real(influencing_cells.size())
  real grey_factor;                 // grey value: 0.: fully liquid, 1.: fully solid
  vec3_t n_vec;                     // approximate normal vector of object boundary
};

struct blackcell_t
{
  size_t l_cell;                // ghost cell index in patch indexing
  vector<size_t> influencing_cells; // cells influencing ghost_cell (on same patch only)
  real average_quote;               // 1./real(influencing_cells.size())
};

class BlockObject
{

  friend class BlockObjectBC;
  friend class CompressibleEulerBOBC;

protected: // data

  PerfectGas m_Gas;
  PatchGrid* m_PatchGrid;
  size_t     m_FieldIndex;

  size_t m_NumBlackLayers;
  size_t m_NumBlackInner;
  size_t m_NumBlackGroups;

  size_t m_GreyResolution;

  ObjectDefinition* m_ObjectDefinition;

  // mem stucture for m_Cells... data;
  //  1st dim: counter index of affected patch
  //  2nd dim: counter index of affected cell in patch
  //  3rd dim (pair::first) : [0]: cell index in patch
  //                          [1]..size() : indices of influencing cells
  //          (pair::second): real (1/(size()-1)
  //OLD  vector<vector<pair<vector<size_t>, real> > > m_CellsFront1;
  //OLD  vector<vector<pair<vector<size_t>, real> > > m_CellsFront2;
  //OLD  vector<vector<pair<vector<size_t>, real> > > m_CellsInside;

  // postponed
  //  1st dim: layer (front0, front, 1, ..., inside) [m_NumBlackLayers + 1]
  //  2nd dim: counter index of affected patch
  //  3rd dim: counter index of affected cell in patch
  //  4th dim (pair::first) : [0]: cell index in patch
  //                          [1]..size() : indices of influencing cells
  //          (pair::second): real (1.(size()-2)
  //NEVER  vector<vector<vector<pair<vector<size_t>, real> > > > m_CellsAllLayers;


  /**
    * Data construction for black cells in a PatchGrid.
    * Black cells are fully inside object(s)
    * mem stucture:
    *   1st dim: layer and inner cell groups:
    *            (front(0), front(1), ... front(m_NumBlackLayers-1), inner(0), inner(1), ... inner(m_NumBlackInner-1))
    *   2nd dim: counter index of affected patch
    *   3rd dim: counter index of affected cell in patch */
  vector<vector<vector<blackcell_t> > > m_BlackCells;

  /**
    * Data construction for grey front cells in a PatchGrid.
    * Grey cells carry individual grey values and are located at the outer front
    * between "black" (solid) and "white" (fluid) cells.
    * mem structure:
    *   1st dim: counter index of affected patch
    *   2nd dim: counter index of affected cell in patch */
  vector<vector<greycell_t> > m_GreyCells;

  /**
    * Indirect indexing of patches affected by this BC.
    * mem structure:
    *   Indexing corresponds to
    *     2nd dimension of m_BlackCells
    *     1st dimension of m_GreyCells */
  vector<size_t> m_AffectedPatchIDs;

  /**
    * Double Indirect indexing of fully black patches pointing onto m_AffectedPatchIDs.
    * These are patches whose core cells are all inside solid object(s). Does not
    * consider eventual grey cells at the patch boundaries, as it is assumed, that at
    * least one seek layer exists (unless for eventual outer boundary conditions).
    * fully black patches may be discarded, if a suitable BC for
    * their remaining neighbour patches exists.
    * Access to a fully black patch:
    *   i_p = m_AffectedPatchIDs[m_FullyBlackPatches[lll]]; */
  vector<size_t> m_FullyBlackPatches;



protected: // methods

  //void copyToHost();
  //void copyToDevice();
  //void processFront(vector<vector<pair<vector<size_t>, real> > >& front, real& p_average, real& T_average);


public: // methods

  BlockObject(PatchGrid* patch_grid,
              size_t num_layers = 2, size_t num_inner = 2,
              size_t grey_resolution = 10);

  void setGreyResolution(real grey_resolution) {m_GreyResolution = grey_resolution;}

  void update(ObjectDefinition* object);

  void analyseGreyCount(Patch* patch, size_t l_cell,
                        size_t& grey_count, size_t& max_grey_count,
                        vec3_t& nxyz);

  void checkRecurrence();

  PatchGrid* getPatchGrid() {return m_PatchGrid; }

  vector<size_t> getAffectedPatchIDs() { return m_AffectedPatchIDs; }

  vector<size_t>* getAffectedPatchIDsPtr() { return &(m_AffectedPatchIDs); }

  vector<vector<vector<blackcell_t> > > getBlackCells() {return m_BlackCells; }

  vector<vector<greycell_t> > getGreyCells() {return m_GreyCells; }

  vector<size_t> getFullyBlackPatches() {return m_FullyBlackPatches; }


  //void setLayerIndexToVar (real** var);
  void setLayerIndexToVar (const size_t &i_field,
                           const size_t &i_variable);

  // must go -> put in BlockBC : public GenericOperation
  //virtual void operator()();
};

#endif // BLOCKOBJECT_H
