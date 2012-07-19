#include "patch.h"

Patch::Patch()
{
  m_Data = NULL;
  m_NumFields = 0;
  m_NumVariables = 0;
  m_VariableSize = 0;
  m_FieldSize = 0;
}

Patch::~Patch()
{
  deleteData();
}

void Patch::allocateData()
{
  m_FieldSize = m_NumVariables * m_VariableSize;
  m_Data = new real [m_NumFields*m_FieldSize];
}

void Patch::deleteData()
{
  delete [] m_Data;
}

void Patch::resize(size_t variable_size)
{
  m_VariableSize = variable_size;
  allocateData();
}

real* Patch::getField(size_t i_field)
{
  return m_Data + i_field*m_FieldSize;
}

real* Patch::getVariable(size_t i_field, size_t i_variable)
{
  return getField(i_field) + i_variable*m_VariableSize;
}

void Patch::computeVariableDifference(size_t i_field1, size_t i_var1, size_t i_field2, size_t i_var2, real &max_norm, real &l2_norm)
{
  restrict real *var1 = getVariable(i_field1, i_var1);
  restrict real *var2 = getVariable(i_field2, i_var2);
  max_norm = 0.0;
  l2_norm  = 0.0;
  for (size_t i = 0; i < m_VariableSize; ++i) {
    real diff = sqr(var2[i] - var1[i]);
    max_norm = max(max_norm, diff);
    l2_norm += diff;
  }
  max_norm = sqrt(max_norm);
  l2_norm  = sqrt(l2_norm);
}
