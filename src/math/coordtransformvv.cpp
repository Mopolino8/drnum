#include "coordtransformvv.h"

CoordTransformVV::CoordTransformVV()
{
  // nothing: m_transform and m_transform_inverse use their defauts
}

void CoordTransformVV::setMatrix(const mat3_t& A)
{
  m_transform.setMatrix(A);
  computeInverse();
}

real CoordTransformVV::setMatrixFromBaseIJ(vec3_t base_i, vec3_t base_j)
{
  real health = m_transform.setMatrixFromBaseIJ(base_i, base_j);
  computeInverse();
  return health;
}

void CoordTransformVV::setVector(const vec3_t& b)
{
  m_transform.setVector(b);
  computeInverse();
}

void CoordTransformVV::setAll(const mat3_t& A, const vec3_t& b)
{
  m_transform.setAll(A, b);
  computeInverse();
}

void CoordTransformVV::setAll(const CoordTransform& transform)
{
  m_transform = transform;
  computeInverse();
}

void CoordTransformVV::scaleVector(real scfactor)
{
  m_transform.scaleVector(scfactor);
  computeInverse();
}

void CoordTransformVV::setTransFromTo(const CoordTransformVV& c_from, const CoordTransformVV& c_to)
{
  m_transform = (c_from.m_transform_inverse).concatenate(c_to.m_transform); // avoids 1st inverting
  m_transform_inverse = (c_to.m_transform_inverse).concatenate(c_from.m_transform); // avoids 2nd inverting
  // m_transform.setTransFromTo(c_from.m_transform, c_to.m_transform);
  // computeInverse();
}

mat3_t CoordTransformVV::getMatrix()
{ 
  mat3_t A = m_transform.getMatrix();
  return A;
}

vec3_t CoordTransformVV::getVector()
{
  return m_transform.getVector();
}

mat3_t CoordTransformVV::getInvMatrix()
{
  mat3_t InvA = m_transform_inverse.getMatrix();;
  return InvA;
}

vec3_t CoordTransformVV::getInvVector()
{
  return m_transform_inverse.getVector();
}

CoordTransform CoordTransformVV::extractForward()
{
  return m_transform;
}

CoordTransform CoordTransformVV::extractReverse()
{
  return m_transform_inverse;
}

void CoordTransformVV::computeInverse()
{
  m_transform_inverse = m_transform.inverse();
}
