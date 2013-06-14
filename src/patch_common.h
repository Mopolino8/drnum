#include "blockcfd.h"

private: // data types

/**
  * Struct holding data for one donor->receiver relation.
  * To be owned by receiving patch.
  */
struct donor_t
{
  size_t variable_size;              ///< number of reals per variable in donor patch
  real* data;                        ///< pointer to m_Data of donor patch (on device being!)
  size_t num_receiver_cells;         ///< Number of cells, receiving data from donor patch
  size_t stride;                     ///< Fixed number of donor cells for each receiving cell
  size_t receiver_index_field_start; ///< Starting index in concatenated receiving cell indicees field of receiving patch
  size_t donor_wi_field_start;       ///< Starting index in concatenated index and weight field for all donor patches
};


private: // attributes

real*   m_Data;         ///< the main data block of this patch
size_t  m_NumFields;    ///< number of fields (e.g. old, new, ...)
size_t  m_NumVariables; ///< number of variables (e.g. rho, rhou, ...)
size_t  m_FieldSize;    ///< length of each field
size_t  m_VariableSize; ///< length of each variable

Transformation m_Transformation; /// @todo merge: kept for compatibility

// Direct transfer lists, omitting InterCoeff classes

size_t m_NumDonorPatches;         ///< Number of donor patches influencing this (receiver) patch
size_t m_NumReceivingCellsConcat; ///< Number of concatenated cells receiving data from any of all donors (multiple indexing)
size_t m_NumReceivingCellsUnique; ///< Number of cells receiving data from any of all donors (unique indexing)
size_t m_NumDonorWIConcat;        ///< Number of concatenated donor cell contributions (all donor patches, all receiving cells times stride)

size_t* m_ReceivingCellIndicesConcat; ///< Concatenated index field of receiving cells in sequence for all donor patches [m_NumReceivingCellsFull]
size_t* m_ReceivingCellIndicesUnique; ///< Index field of receiving cells in unique sequence [m_NumReceivingCellsUnique]
/// @todo m_ReceivingCellIndexUnique equivalent to m_receive_cells. Clean up!!!

donor_t* m_Donors;         ///< All donor data structs for this (receiver) patch

size_t* m_DonorIndexConcat;  ///< Concatenated donor cell indicees [m_NumDonorWIConcat]
real*  m_DonorWeightConcat;  ///< Concatenated donor cell weights [m_NumDonorWIConcat]
// that would be nicer: pair<size_t, real>*  m_DonorWI;


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



