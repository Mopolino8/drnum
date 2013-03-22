#include "blockcfd.h"

private:

real*   m_Data;         ///< the main data block of this patch
size_t  m_NumFields;    ///< number of fields (e.g. old, new, ...)
size_t  m_NumVariables; ///< number of variables (e.g. rho, rhou, ...)
size_t  m_FieldSize;    ///< length of each field
size_t  m_VariableSize; ///< length of each variable

Transformation m_Transformation; /// @todo merge: kept for compatibility


public:

CUDA_DH real* getField(size_t i_field)
{
  return m_Data + i_field*m_FieldSize;
}

CUDA_DH real* getVariable(size_t i_field, size_t i_variable)
{
  return getField(i_field) + i_variable*m_VariableSize;
}

/**
  * This is one of the main methods which will be used for data exchange.
  * How this can be used to handle exchange between CPU/GPU, CPU/CPU, GPU/GPU, and NODE/NODE (MPI) needs
  * to be established as soon as possible.
  * @todo look into different data exchange methods
  * @param i_field the index of the field (e.g. old, new, ...)
  * @param i_var the index of the variable (e.g. rho, rhou, ...)
  * @param i the spatial index (either cell or node)
  * @return the value at the specified position
  */
CUDA_DH real getValue(size_t i_field, size_t i_var, size_t i)
{
  return m_Data[i_field*m_FieldSize + i_var*m_VariableSize + i];
}

CUDA_DH size_t numFields() const
{
  return m_NumFields;
}

CUDA_DH size_t numVariables() const
{
  return m_NumVariables;
}

CUDA_DH size_t dataSize() const
{
  return m_NumFields*m_FieldSize;
}

CUDA_DH size_t fieldSize() const
{
  return m_FieldSize;
}

CUDA_DH size_t variableSize() const
{
  return m_VariableSize;
}

CUDA_DH real* getData()
{
  return m_Data;
}

/**
 * Copy simple data attributes from another object.
 * The other object can have a different type as long as the required attributes are present.
 * param obj a constant reference to the other object
 */
template <typename T>
CUDA_HO void copyAttributes(T* obj)
{
  m_NumFields    = obj->numFields();
  m_NumVariables = obj->numVariables();
  m_FieldSize    = obj->fieldSize();
  m_VariableSize = obj->variableSize();
}



