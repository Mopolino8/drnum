#ifndef PATCH_H
#define PATCH_H

#include <cstddef>

#include "blockcfd.h"

class Patch
{

protected: // attributes

  real*   m_Data;      ///< the main data block of this patch
  size_t  m_NumFields; ///< number of fields (e.g. rho, rhou, ...)
  size_t  m_FieldSize; ///< length of each field

protected: // methods

  void  allocateData(size_t field_size);
  real* getField(size_t i_field);

public: // methods

  Patch();
  virtual ~Patch();

  /**
    * This is one of the main methods which will be used for data exchange.
    * How this can be used to handle exchange between CPU/GPU, CPU/CPU, GPU/GPU, and NODE/NODE (MPI) needs
    * to be established as sooon as possible.
    * @todo look into different data exchange methods
    * @param i_field the index of the field (e.g. rho, rhou, ...)
    * @param i the spatial index (either cell or node)
    * @return the value at the specified position
    */
  real getValue(size_t i_field, size_t i) { return m_Data[i_field*m_FieldSize + i]; }

};

#endif // PATCH_H
