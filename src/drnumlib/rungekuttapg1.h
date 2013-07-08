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
