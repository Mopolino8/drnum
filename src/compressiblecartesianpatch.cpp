#include "compressiblecartesianpatch.h"

CompressibleCartesianPatch::CompressibleCartesianPatch()
{
  setNumberOfFields(15);
}

void CompressibleCartesianPatch::updatePointers()
{
  m_Rho     = getField(0);
  m_Rhou    = getField(1);
  m_Rhov    = getField(2);
  m_Rhow    = getField(3);
  m_RhoE    = getField(4);
  m_OldRho  = getField(5);
  m_OldRhou = getField(6);
  m_OldRhov = getField(7);
  m_OldRhow = getField(8);
  m_OldRhoE = getField(9);
  m_ResRho  = getField(10);
  m_ResRhou = getField(11);
  m_ResRhov = getField(12);
  m_ResRhow = getField(13);
  m_ResRhoE = getField(14);
}

void CompressibleCartesianPatch::computeAusmPlus()
{
  restrict real* r  = rho();
  restrict real* ru = rhou();
  restrict real* rv = rhov();
  restrict real* rw = rhow();
  restrict real* rE = rhoE();
  restrict real* res_r  = resRho();
  restrict real* res_ru = resRhou();
  restrict real* res_rv = resRhov();
  restrict real* res_rw = resRhow();
  restrict real* res_rE = resRhoE();
}
