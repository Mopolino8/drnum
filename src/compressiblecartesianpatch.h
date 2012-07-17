#ifndef COMPRESSIBLECARTESIANPATCH_H
#define COMPRESSIBLECARTESIANPATCH_H

#include "cartesianpatch.h"

class CompressibleCartesianPatch : public CartesianPatch
{

private: // attributes

  real* m_OldRho;
  real* m_OldRhou;
  real* m_OldRhov;
  real* m_OldRhow;
  real* m_OldRhoE;
  real* m_Rho;
  real* m_Rhou;
  real* m_Rhov;
  real* m_Rhow;
  real* m_RhoE;
  real* m_ResRho;
  real* m_ResRhou;
  real* m_ResRhov;
  real* m_ResRhow;
  real* m_ResRhoE;

protected: // methods

  virtual void updatePointers();

  real* oldRho  () { return m_OldRho; }
  real* oldRhou () { return m_OldRhou; }
  real* oldRhov () { return m_OldRhov; }
  real* oldRhow () { return m_OldRhow; }
  real* oldRhoE () { return m_OldRhoE; }
  real* rho     () { return m_Rho; }
  real* rhou    () { return m_Rhou; }
  real* rhov    () { return m_Rhov; }
  real* rhow    () { return m_Rhow; }
  real* rhoE    () { return m_RhoE; }
  real* resRho  () { return m_ResRho; }
  real* resRhou () { return m_ResRhou; }
  real* resRhov () { return m_ResRhov; }
  real* resRhow () { return m_ResRhow; }
  real* resRhoE () { return m_ResRhoE; }

public: // methods

  CompressibleCartesianPatch();
  void computeAusmPlus();

};

#endif // COMPRESSIBLECARTESIANPATCH_H
