#include "cartesianpatch.h"

CartesianPatch::CartesianPatch()
{
  m_NumI = 1;
  m_NumJ = 1;
  m_NumK = 1;
}

void CartesianPatch::computeDeltas()
{
  double Lx = sqrt(sqr(m_Uxo) + sqr(m_Uyo) + sqr(m_Uzo));
  double Ly = sqrt(sqr(m_Vxo) + sqr(m_Vyo) + sqr(m_Vzo));
  double Lz = sqrt(sqr(m_Wxo) + sqr(m_Wyo) + sqr(m_Wzo));
  countFlops(15);
  countSqrts(3);

  m_Dx = Lx/m_NumI;
  m_Dy = Ly/m_NumJ;
  m_Dz = Lz/m_NumK;
  countFlops(3);

  m_InvDx = 1.0/m_Dx;
  m_InvDy = 1.0/m_Dx;
  m_InvDz = 1.0/m_Dx;
  countFlops(3);

  m_Dixo = m_Uxo/m_NumI;
  m_Diyo = m_Uyo/m_NumI;
  m_Dizo = m_Uzo/m_NumI;
  countFlops(3);

  m_Djxo = m_Vxo/m_NumJ;
  m_Djyo = m_Vyo/m_NumJ;
  m_Djzo = m_Vzo/m_NumJ;
  countFlops(3);

  m_Djxo = m_Wxo/m_NumK;
  m_Djyo = m_Wyo/m_NumK;
  m_Djzo = m_Wzo/m_NumK;
  countFlops(3);

  m_Xco = m_Xo + 0.5*(m_Dixo + m_Djxo + m_Dkxo);
  m_Yco = m_Yo + 0.5*(m_Diyo + m_Djyo + m_Dkyo);
  m_Zco = m_Zo + 0.5*(m_Dizo + m_Djzo + m_Dkzo);
  countFlops(12);
}

void CartesianPatch::setupAligned(real x1, real y1, real z1, real x2, real y2, real z2)
{
  m_Xo = x1;
  m_Yo = y1;
  m_Zo = z1;

  m_Uxo = x2 - x1;
  m_Uyo = 0;
  m_Uzo = 0;
  countFlops(1);

  m_Vxo = 0;
  m_Vyo = y2 - y1;
  m_Vzo = 0;
  countFlops(1);

  m_Wxo = 0;
  m_Wyo = 0;
  m_Wzo = z2 - z1;
  countFlops(1);

  computeDeltas();

  Transformation t;
  t.setVector(vec3_t(x1, y1, z1));
  setTransformation(t.inverse());
}

void CartesianPatch::resize(size_t num_i, size_t num_j, size_t num_k)
{
  m_NumI = num_i;
  m_NumJ = num_j;
  m_NumK = num_k;
  deleteData();
  Patch::resize(m_NumI*m_NumJ*m_NumK);
  computeDeltas();
}

