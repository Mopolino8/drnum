#include "patch.h"

Patch::Patch()
{
  m_Data = NULL;
  m_NumFields = 0;
}

Patch::~Patch()
{
  delete [] m_Data;
}

void Patch::allocateData(size_t field_size)
{
  m_FieldSize = field_size;
  m_Data = new real [m_NumFields*m_FieldSize];
}

real* Patch::getField(size_t i_field)
{
  return m_Data + i_field*m_FieldSize;
}
