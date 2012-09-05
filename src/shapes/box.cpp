#include "box.h"

Box::Box()
{
}

void Box::setGeometry(real x1, real y1, real z1, real x2, real y2, real z2)
{
  m_X1 = x1;
  m_Y1 = y1;
  m_Z1 = z1;
  m_X2 = x2;
  m_Y2 = y2;
  m_Z2 = z2;
}

void Box::transform(const Transformation &transformation)
{
  m_Transform = transformation;
  m_InvTransform = m_Transform.inverse();
}

void Box::reset()
{
  transform(Transformation());
}

bool Box::getBoundaryMetric(real x1, real y1, real z1,
                            real x2, real y2, real z2,
                            real &k,
                            real &nx, real &ny, real &nz)
{
  /// @todo implement rotation of normal vectors

  m_InvTransform(x1, y1, z1);
  m_InvTransform(x2, y2, z2);

  // left
  if ((x1 <= m_X1 && x2 > m_X1) || (x1 >= m_X1 && x2 < m_X1)) {
    k = (m_X1 - x1)/(x2 - x1);
    real y = y1 + k*(y2 - y1);
    real z = z1 + k*(z2 - z1);
    if (y > m_Y1 && y < m_Y2 && z > m_Z1 && z < m_Z2) {
      nx = -1;
      ny = 0;
      nz = 0;
      return true;
    }
  }

  // right
  if ((x1 <= m_X1 && x2 > m_X1) || (x1 >= m_X1 && x2 < m_X1)) {
    k = (m_X2 - x2)/(x1 - x2);
    real y = y1 + k*(y2 - y1);
    real z = z1 + k*(z2 - z1);
    if (y > m_Y1 && y < m_Y2 && z > m_Z1 && z < m_Z2) {
      nx = 1;
      ny = 0;
      nz = 0;
      return true;
    }
  }

  // front
  if ((y1 <= m_Y1 && y2 > m_Y1) || (y1 >= m_Y1 && y2 < m_Y1)) {
    k = (m_Y1 - y1)/(y2 - y1);
    real x = x1 + k*(x2 - x1);
    real z = z1 + k*(z2 - z1);
    if (x > m_X1 && x < m_X2 && z > m_Z1 && z < m_Z2) {
      nx = 0;
      ny = -1;
      nz = 0;
      return true;
    }
  }

  // back
  if ((y1 <= m_Y1 && y2 > m_Y1) || (y1 >= m_Y1 && y2 < m_Y1)) {
    k = (m_Y2 - y2)/(y1 - y2);
    real x = x1 + k*(x2 - x1);
    real z = z1 + k*(z2 - z1);
    if (x > m_X1 && x < m_X2 && z > m_Z1 && z < m_Z2) {
      nx = 0;
      ny = 1;
      nz = 0;
      return true;
    }
  }

  // bottom
  if ((z1 <= m_Z1 && z2 > m_Z1) || (z1 >= m_Z1 && x2 < m_Z1)) {
    k = (m_Z1 - z1)/(z2 - z1);
    real x = x1 + k*(x2 - x1);
    real y = y1 + k*(y2 - y1);
    if (x > m_X1 && x < m_X2 && y > m_Y1 && y < m_Y2) {
      nx = 0;
      ny = 0;
      nz = -1;
      return true;
    }
  }

  // top
  if ((z1 <= m_Z2 && z2 > m_Z2) || (z1 >= m_Z2 && z2 < m_Z2)) {
    k = (m_Z2 - z2)/(z1 - z2);
    real x = x1 + k*(x2 - x1);
    real y = y1 + k*(y2 - y1);
    if (x > m_X1 && x < m_X2 && y > m_Y1 && y < m_Y2) {
      nx = 0;
      ny = 0;
      nz = 1;
      return true;
    }
  }

  return false;
}

bool Box::isInside(real x, real y, real z)
{
  m_InvTransform(x, y, z);
  if (x > m_X1 && x < m_X2 && y > m_Y1 && y < m_Y2 && z > m_Z1 && z < m_Z2) {
    return true;
  }
  return false;
}
