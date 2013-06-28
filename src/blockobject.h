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

  // 1st dim: patch id, 2nd dim: cell in patch , 3rd dim: neighbour cells to take data from
  vector<vector<pair<vector<size_t>, real> > > m_CellsFront1;
  vector<vector<pair<vector<size_t>, real> > > m_CellsFront2;
  vector<vector<pair<vector<size_t>, real> > > m_CellsInside;

  vector<size_t> m_affectedPatchIDs;

public: // methods

  BlockObject();

  void fillAll();

  virtual bool isInside(const real& xo, const real& yo, const real& zo) = 0;

  bool isInside(vec3_t xyzo);

  vector<size_t> getAffectedPatchIDs() {
    return m_affectedPatchIDs;
  }

  vector<size_t>* getAffectedPatchIDsPtr() {
    return &(m_affectedPatchIDs);
  }

  virtual void operator ()();

  void applyFrontBC (const Patch* patch, const size_t& l);

  void applyInnerBC (const Patch* patch, const size_t& l);

};

inline bool BlockObject::isInside(vec3_t xyzo)
{
  return isInside (xyzo[0], xyzo[1],xyzo[2]);
}

#endif // BLOCKOBJECT_H
