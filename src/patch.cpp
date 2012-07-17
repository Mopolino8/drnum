#include "patch.h"

Patch::Patch()
{
  m_Data = NULL;
  m_NumFields = 0;
}

Patch::~Patch()
{
  deleteData();
}

void Patch::allocateData(size_t field_size)
{
  m_FieldSize = field_size;
  m_Data = new real [m_NumFields*m_FieldSize];
  updatePointers();
}

void Patch::deleteData()
{
  delete [] m_Data;
}

real* Patch::getField(size_t i_field)
{
  return m_Data + i_field*m_FieldSize;
}

void Patch::computeFieldDifference(size_t i_field1, size_t i_field2, real &max_norm, real &l2_norm)
{
  restrict real *field1 = getField(i_field1);
  restrict real *field2 = getField(i_field2);
  max_norm = 0.0;
  l2_norm = 0.0;
  for (size_t i = 0; i < m_FieldSize; ++i) {
    real diff = sqr(field2[i] - field1[i]);
    max_norm = max(max_norm, diff);
    l2_norm += diff;
  }
  max_norm = sqrt(max_norm);
  l2_norm = sqrt(l2_norm);
}
