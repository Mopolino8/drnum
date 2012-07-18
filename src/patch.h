#ifndef PATCH_H
#define PATCH_H

#include <cstddef>

#include "blockcfd.h"

class Patch
{

private: // attributes

  real*   m_Data;      ///< the main data block of this patch
  size_t  m_NumFields; ///< number of fields (e.g. rho, rhou, ...)
  size_t  m_FieldSize; ///< length of each field

protected: // methods

  void  allocateData(size_t field_size);
  void  deleteData();
  real* getField(size_t i_field);
  void  setNumberOfFields(size_t num_fields) { m_NumFields = num_fields; }
  void  setFieldToZero(real *field);
  void  copyField(real *src, real *dst);

  /**
   * Will be called at the end of allocaData().
   * This can be used to update pointers after a resize operation.
   */
  virtual void updatePointers() {}


public: // methods

  Patch();
  virtual ~Patch();

  /**
    * This is one of the main methods which will be used for data exchange.
    * How this can be used to handle exchange between CPU/GPU, CPU/CPU, GPU/GPU, and NODE/NODE (MPI) needs
    * to be established as soon as possible.
    * @todo look into different data exchange methods
    * @param i_field the index of the field (e.g. rho, rhou, ...)
    * @param i the spatial index (either cell or node)
    * @return the value at the specified position
    */
  real getValue(size_t i_field, size_t i) { return m_Data[i_field*m_FieldSize + i]; }

  /**
   * Compute the difference between two fields. This is intended to be used for convergence monitoring.
   * @param i_field1 index of the first field
   * @param i_field2 index of the first field
   * @param max_norm will hold the maximal absolute difference
   * @param l2_norm will hold the L2 norm of the difference
   */
  void computeFieldDifference(size_t i_field1, size_t i_field2, real &max_norm, real &l2_norm);

  size_t numFields() { return m_NumFields; }
  size_t fieldSize() { return m_FieldSize; }

  virtual void computeResiduals() = 0;
  virtual void subStep(real dt) = 0;

};

inline void Patch::setFieldToZero(real *field)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    field[i] = 0.0;
  }
}

inline void Patch::copyField(real *src, real *dst)
{
  for (size_t i = 0; i < m_FieldSize; ++i) {
    dst[i] = src[i];
  }
}

#endif // PATCH_H
