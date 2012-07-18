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

void CompressibleCartesianPatch::setState(size_t i, size_t j, size_t k, real p, real T)
{
  #include "code_blocks/compressible_variables.h"

  setF(r,  i, j, k, p/(gasR() * T));
  setF(ru, i, j, k, 0);
  setF(rv, i, j, k, 0);
  setF(rw, i, j, k, 0);
  setF(rE, i, j, k, p/(gasGamma() - 1));
}

void CompressibleCartesianPatch::setOldState()
{
  copyField(m_Rho,  m_OldRho);
  copyField(m_Rhou, m_OldRhou);
  copyField(m_Rhov, m_OldRhov);
  copyField(m_Rhow, m_OldRhow);
  copyField(m_RhoE, m_OldRhoE);
}

void CompressibleCartesianPatch::subStep(real dt)
{
  #include "code_blocks/compressible_variables.h"
  #include "code_blocks/compressible_residuals.h"

  setFieldToZero(res_r);
  setFieldToZero(res_ru);
  setFieldToZero(res_rv);
  setFieldToZero(res_rw);
  setFieldToZero(res_rE);
  computeResiduals();
  for (size_t i = 0; i < sizeI(); ++i) {
    for (size_t j = 0; j < sizeJ(); ++j) {
      for (size_t k = 0; k < sizeK(); ++k) {
        setF(r,  i, j, k, f(m_OldRho,  i, j, k) + dt*f(res_r,  i, j, k));
        setF(ru, i, j, k, f(m_OldRhou, i, j, k) + dt*f(res_ru, i, j, k));
        setF(rv, i, j, k, f(m_OldRhov, i, j, k) + dt*f(res_rv, i, j, k));
        setF(rw, i, j, k, f(m_OldRhow, i, j, k) + dt*f(res_rw, i, j, k));
        setF(rE, i, j, k, f(m_OldRhoE, i, j, k) + dt*f(res_rE, i, j, k));
      }
    }
  }
}
