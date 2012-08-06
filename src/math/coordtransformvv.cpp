#include "coordtransformvv.h"

CoordTransformVV::CoordTransformVV()
{
  // nothing: m_transform and m_transform_inverse use their defauts
}

void CoordTransformVV::setMatrix(mat3_t A)
{
  m_transform.setMatrix(A);
  computeInverse();
}

void CoordTransformVV::setVector(vec3_t b)
{
  m_transform.setVector(b);
  computeInverse();
}

void CoordTransformVV::setAll(mat3_t A, vec3_t b)
{
  m_transform.setAll(A, b);
  computeInverse();
}

void CoordTransformVV::setTransFromTo(const CoordTransformVV& c_from, const CoordTransformVV& c_to)
{
  m_transform.setTransFromTo(c_from.m_transform, c_to.m_transform);
  computeInverse();
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

void CoordTransformVV::computeInverse()
{
  m_transform_inverse = m_transform.inverse();
}
