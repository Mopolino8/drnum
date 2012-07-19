#include "compressibleflux.h"

CompressibleFlux::CompressibleFlux()
{
  m_Rho  = NULL;
  m_Rhou = NULL;
  m_Rhov = NULL;
  m_Rhow = NULL;
  m_RhoE = NULL;
  m_ResRho  = NULL;
  m_ResRhou = NULL;
  m_ResRhov = NULL;
  m_ResRhow = NULL;
  m_ResRhoE = NULL;
}

