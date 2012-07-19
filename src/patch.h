#ifndef PATCH_H
#define PATCH_H

#include <cstddef>

#include "blockcfd.h"

class Patch
{

private: // attributes

  real*   m_Data;         ///< the main data block of this patch
  size_t  m_NumFields;    ///< number of fields (e.g. old, new, ...)
  size_t  m_NumVariables; ///< number of variables (e.g. rho, rhou, ...)
  size_t  m_FieldSize;    ///< length of each field
  size_t  m_VariableSize; ///< length of each variable


private: // methods

  void  allocateData();


protected: // methods

  void  deleteData();
  void  resize(size_t variable_size);

  real* getField(size_t i_field);
  real* getVariable(size_t i_field, size_t i_variable);
  void  setNumberOfFields(size_t num_fields) { m_NumFields = num_fields; }
  void  setNumberOfVariables(size_t num_variables) { m_NumVariables = num_variables; }


public: // methods

  Patch();
  virtual ~Patch();

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
  real getValue(size_t i_field, size_t i_var, size_t i) { return m_Data[i_field*m_FieldSize + i_var*m_VariableSize + i]; }

  /**
   * Compute the difference between two variables. This is intended to be used for convergence monitoring.
   * @param i_field1 index of the first field
   * @param i_var1 index of the first variable
   * @param i_field2 index of the first field
   * @param i_var2 index of the second variable
   * @param max_norm will hold the maximal absolute difference
   * @param l2_norm will hold the L2 norm of the difference
   */
  void computeVariableDifference(size_t i_field1, size_t i_var1, size_t i_field2, size_t i_var2, real &max_norm, real &l2_norm);

  size_t numFields()    { return m_NumFields; }
  size_t numVariables() { return m_NumVariables; }
  size_t fieldSize()    { return m_FieldSize; }
  size_t variableSize() { return m_VariableSize; }

  virtual void preStep() {}
  virtual void postStep() {}
  virtual void subStep(real dt) = 0;

  void  setFieldToZero(real *field);
  void  setFieldToZero(size_t i_field);
  void  copyField(real *src, real *dst);
  void  copyField(size_t i_src, size_t i_dst);
  void  addField(real *src, real factor, real *dst);
  void  addField(real *op1, real factor, real *op2, real *dst);
  void  addField(size_t i_src, size_t factor, size_t i_dst);
  void  addField(size_t i_op1, size_t factor, size_t i_op2, size_t i_dst);

};

inline void Patch::setFieldToZero(real *field)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    field[i] = 0.0;
  }
}

inline void Patch::setFieldToZero(size_t i_field)
{
  setFieldToZero(getField(i_field));
}

inline void Patch::copyField(real *src, real *dst)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    dst[i] = src[i];
  }
}

inline void Patch::copyField(size_t i_src, size_t i_dst)
{
  copyField(getField(i_src), getField(i_dst));
}

inline void Patch::addField(real *src, real factor, real *dst)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    dst[i] += factor*src[i];
  }
}

inline void Patch::addField(real *op1, real factor, real *op2, real *dst)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    dst[i] = op1[i] + factor*op2[i];
  }
}

inline void Patch::addField(size_t i_src, size_t factor, size_t i_dst)
{
  addField(getField(i_src), factor, getField(i_dst));
}

inline void Patch::addField(size_t i_op1, size_t factor, size_t i_op2, size_t i_dst)
{
  addField(getField(i_op1), factor, getField(i_op2), getField(i_dst));
}

#endif // PATCH_H
