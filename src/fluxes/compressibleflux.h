#ifndef COMPRESSIBLEFLUX_H
#define COMPRESSIBLEFLUX_H

#include "compressiblecartesianpatch.h"

class CompressibleFlux : public CompressibleObject
{

protected: // attributes

  restrict real *m_Rho;
  restrict real *m_Rhou;
  restrict real *m_Rhov;
  restrict real *m_Rhow;
  restrict real *m_RhoE;
  restrict real *m_ResRho;
  restrict real *m_ResRhou;
  restrict real *m_ResRhov;
  restrict real *m_ResRhow;
  restrict real *m_ResRhoE;

public:

  CompressibleFlux();

  void setRho (real *r)  { m_Rho  = r; }
  void setRhou(real *ru) { m_Rhou = ru; }
  void setRhov(real *rv) { m_Rhov = rv; }
  void setRhow(real *rw) { m_Rhow = rw; }
  void setRhoE(real *rE) { m_RhoE = rE; }

  void setResRho (real *res_r)  { m_ResRho  = res_r; }
  void setResRhou(real *res_ru) { m_ResRhou = res_ru; }
  void setResRhov(real *res_rv) { m_ResRhov = res_rv; }
  void setResRhow(real *res_rw) { m_ResRhow = res_rw; }
  void setResRhoE(real *res_rE) { m_ResRhoE = res_rE; }

};

#endif // COMPRESSIBLEFLUX_H
