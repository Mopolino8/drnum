#include "shapes/halfspace.h"

void HalfSpace::transform(const Transformation &transformation)
{
  // @todo look at generality (e.g. this would fall over for a rotation)
  m_XOrig = m_X;
  m_YOrig = m_Y;
  m_ZOrig = m_Z;
  m_NxOrig = m_Nx;
  m_NyOrig = m_Ny;
  m_NzOrig = m_Nz;

  transformation(m_X, m_Y, m_Z);
}

void HalfSpace::reset()
{
  m_X = m_XOrig;
  m_Y = m_YOrig;
  m_Z = m_ZOrig;
  m_Nx = m_NxOrig;
  m_Ny = m_NyOrig;
  m_Nz = m_NzOrig;
}
