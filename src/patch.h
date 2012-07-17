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

};

#endif // PATCH_H
