// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#include "drnum.h"

private: // attributes

real*   m_Data;         ///< the main data block of this patch
size_t  m_NumFields;    ///< number of fields (e.g. old, new, ...)
size_t  m_NumVariables; ///< number of variables (e.g. rho, rhou, ...)
size_t  m_FieldSize;    ///< length of each field
size_t  m_VariableSize; ///< length of each variable

PatchGrid* m_PatchGrid; ///< the patch grid this patch belongs to

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


// data for split faces (immersed boundary method)
size_t       m_NumSplitFaces;
bool        *m_IsInsideCell;
bool        *m_IsSplitCell;
splitface_t *m_SplitFaces;

protected: // attributes

size_t  m_MyIndex;      ///< Index of patch in sequence of PatchGrid::m_patches. Optional setting.


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

CUDA_DH size_t getNumDonorPatches()
{
  return m_NumDonorPatches;
}

CUDA_DH size_t getNumReceivingCellsConcat()
{
  return m_NumReceivingCellsConcat;
}

CUDA_DH size_t getNumReceivingCellsUnique()
{
  return m_NumReceivingCellsUnique;
}

CUDA_DH size_t getNumDonorWIConcat()
{
  return m_NumDonorWIConcat;
}

CUDA_DH size_t* getReceivingCellIndicesConcat()
{
  return m_ReceivingCellIndicesConcat;
}

CUDA_DH size_t* getReceivingCellIndicesUnique()
{
  return m_ReceivingCellIndicesUnique;
}

CUDA_DH donor_t* getDonors()
{
  return m_Donors;
}

CUDA_DH size_t* getDonorIndexConcat()
{
  return m_DonorIndexConcat;
}

CUDA_DH real* getDonorWeightConcat()
{
  return m_DonorWeightConcat;
}

CUDA_HO PatchGrid* getPatchGrid()
{
  return m_PatchGrid;
}

/**
 * @brief Set a variable set at a specified l_c - position.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param varset the variable set
 */
CUDA_DH void setVarset(size_t i_field, size_t l_c, real* varset)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    *(start_p + i_var*m_VariableSize) = varset[i_var];
    //    m_Data[i_field*m_FieldSize + i_var*m_VariableSize + l_c] = varset[i_var];
  }
}


/**
 * @brief Set a variable set at a specified l_c - position.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param varset the variable set
 * @param start_iv starting var-index (example: 0 for rho, etc...)
 * @param after_last_iv after last var-index (example: 2 for rho_v, etc...)
 * @param start_ivss starting subset var-index (e.g. 0 for var[0])
 */
CUDA_DH void setVarsubset(size_t i_field, size_t l_c, real* varset,
                          size_t start_iv, size_t after_last_iv,
                          size_t start_ivss = 0)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = start_iv; i_var < after_last_iv; ++i_var) {
    size_t i_var_ss = i_var - start_iv + start_ivss;
    *(start_p + i_var*m_VariableSize) = varset[i_var_ss];
  }
}


/**
 * @brief Set a variable set to at a specified l_c - position to (0,0,0...)
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 */
CUDA_DH void setVarsetToZero(size_t i_field, size_t l_c)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    *(start_p + i_var*m_VariableSize) = 0.;
    //    m_Data[i_field*m_FieldSize + i_var*m_VariableSize + l_c] = 0.;
  }
}


/**
 * @brief Set a variable set to at a specified l_c - position to (0,0,0...)
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param start_iv starting var-index (example: 0 for rho, etc...)
 * @param after_last_iv after last var-index (example: 2 for rho_v, etc...)
 */
CUDA_DH void setVarsubsetToZero(size_t i_field, size_t l_c,
                                size_t start_iv, size_t after_last_iv)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = start_iv; i_var < after_last_iv; ++i_var) {
    *(start_p + i_var*m_VariableSize) = 0.;
    //    m_Data[i_field*m_FieldSize + i_var*m_VariableSize + l_c] = 0.;
  }
}


/**
 * @brief Add a varset to a specified varset at an l_c - position.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param varset the variable set
 */
CUDA_DH void addToVarset(size_t i_field, size_t l_c, real* varset)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    *(start_p + i_var*m_VariableSize) += varset[i_var];
  }
}


/**
 * @brief Add a varset to a specified varset at an l_c - position.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param varset the variable set
 * @param start_iv starting var-index (example: 0 for rho, etc...)
 * @param after_last_iv after last var-index (example: 2 for rho_v, etc...)
 * @param start_ivss starting subset var-index (e.g. 0 for var[0])
 */
CUDA_DH void addToVarsubset(size_t i_field, size_t l_c, real* varset,
                            size_t start_iv, size_t after_last_iv,
                            size_t start_ivss = 0)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = start_iv; i_var < after_last_iv; ++i_var) {
    size_t i_var_ss = i_var - start_iv + start_ivss;
    *(start_p + i_var*m_VariableSize) += varset[i_var_ss];
  }
}


/**
 * @brief Add a varset at position l_give to varset at position l_rec .
 * @param i_field the field index
 * @param l_rec the receiving cell index (in 1D-alignment, e.g. unstructured)
 * @param l_give the "giving" cell index.
 */
CUDA_DH void addToVarset(size_t i_field, size_t l_rec, size_t l_give)
{
  real* start_rec = m_Data + i_field*m_FieldSize + l_rec; // for i_var = 0;
  real* start_give = m_Data + i_field*m_FieldSize + l_give; // for i_var = 0;
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    *(start_rec + i_var*m_VariableSize) += *(start_give + i_var*m_VariableSize);
  }
}


/**
 * @brief Add a varset at position l_give to varset at position l_rec .
 * @param i_field the field index
 * @param l_rec the receiving cell index (in 1D-alignment, e.g. unstructured)
 * @param l_give the "giving" cell index.
 * @param start_iv starting var-index (example: 0 for rho, etc...)
 * @param after_last_iv after last var-index (example: 2 for rho_v, etc...)
 */
CUDA_DH void addToVarsubset(size_t i_field, size_t l_rec, size_t l_give,
                            size_t start_iv, size_t after_last_iv)
{
  real* start_rec = m_Data + i_field*m_FieldSize + l_rec; // for i_var = 0;
  real* start_give = m_Data + i_field*m_FieldSize + l_give; // for i_var = 0;
  for (size_t i_var = start_iv; i_var < after_last_iv; ++i_var) {
    *(start_rec + i_var*m_VariableSize) += *(start_give + i_var*m_VariableSize);
  }
}


/**
 * @brief Multiply a varset at position l_c with a scalar value.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param scalar the scalar value to multiply with
 */
CUDA_DH void multVarsetScalar(const size_t& i_field, const size_t& l_c, const real& scalar)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    *(start_p + i_var*m_VariableSize) *= scalar;
  }
}


/**
 * @brief Multiply a varset at position l_c with a scalar value.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param scalar the scalar value to multiply with
 * @param start_iv starting var-index (example: 0 for rho, etc...)
 * @param after_last_iv after last var-index (example: 2 for rho_v, etc...)
 */
CUDA_DH void multVarsubsetScalar(const size_t& i_field, const size_t& l_c, const real& scalar,
                                 size_t start_iv, size_t after_last_iv)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = start_iv; i_var < after_last_iv; ++i_var) {
    *(start_p + i_var*m_VariableSize) *= scalar;
  }
}


/**
 * @brief Get a varset from l_c - position.
 *        Better avoid this function due to its mem jumping.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param ret_varset the conservative variable set (return: write after pointer)
 */
CUDA_DH real* getVarset(size_t i_field, size_t l_c, real* ret_varset)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = 0; i_var < numVariables(); ++i_var) {
    ret_varset[i_var] = *(start_p + i_var*m_VariableSize);
  }
}


/**
 * @brief Get a varset from l_c - position.
 *        Better avoid this function due to its mem jumping.
 * @param i_field the field index
 * @param l_c the cell index (in 1D-alignment, e.g. unstructured)
 * @param ret_varset the conservative variable set (return: write after pointer)
 * @param start_iv starting var-index (example: 0 for rho, etc...)
 * @param after_last_iv after last var-index (example: 2 for rho_v, etc...)
 * @param start_ivss starting subset var-index (e.g. 0 for var[0])
 */
CUDA_DH real* getVarsubset(size_t i_field, size_t l_c, real* ret_varset,
                           size_t start_iv, size_t after_last_iv,
                           size_t start_ivss = 0)
{
  real* start_p = m_Data + i_field*m_FieldSize + l_c; // for i_var = 0;
  for (size_t i_var = start_iv; i_var < after_last_iv; ++i_var) {
    size_t i_var_ss = i_var - start_iv + start_ivss;
    ret_varset[i_var_ss] = *(start_p + i_var*m_VariableSize);
  }
}

/**
 * @brief Get a variable set at a specified index.
 * @param i_field the field index
 * @param i the cell index
 * @param var will hold the conservative variable set afterwards (needs to be allocated beforehand)
 */
template <typename DIM>
CUDA_DH void getVar(DIM, size_t i_field, size_t i, real *var)
{
  for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
    var[i_var] = getVariable(i_field, i_var)[i];
  }
}

/**
 * @brief Get a variable set at a specified index. This method take the dimension (number of variables)
 *        as a real parameter and should not be used in performance critical routines.
 * @param dim the number of variables
 * @param i_field the field index
 * @param i the cell index
 * @param var will hold the conservative variable set afterwards (needs to be allocated beforehand)
 */
CUDA_DH void getVarDim(unsigned int dim, size_t i_field, size_t i, real *var)
{
  for (size_t i_var = 0; i_var < dim; ++i_var) {
    var[i_var] = getVariable(i_field, i_var)[i];
  }
}

/**
 * @brief Set a variable set at a specified index.
 * @param i_field the field index
 * @param i the cell index
 * @param the conservative variable set
 */
template <typename DIM>
CUDA_DH void setVar(DIM, size_t i_field, size_t i, real* var)
{
  for (size_t i_var = 0; i_var < DIM::dim; ++i_var) {
    getVariable(i_field, i_var)[i] = var[i_var];
  }
}

/**
 * @brief Set a variable set at a specified index. This method take the dimension (number of variables)
 *        as a real parameter and should not be used in performance critical routines.
 * @param dim the number of variables
 * @param i_field the field index
 * @param i the cell index
 * @param the conservative variable set
 */
CUDA_DH void setVarDim(unsigned int dim, size_t i_field, size_t i, real *var)
{
  for (size_t i_var = 0; i_var < dim; ++i_var) {
    getVariable(i_field, i_var)[i] = var[i_var];
  }
}

/**
  * Access
  * @return the index of the patch in sequence of PatchGrid::m_patches
  */
CUDA_DH size_t getIndex() {return m_MyIndex;}


/**
 * @brief Get the total number of split faces for the immersed boundaries.
 * @return number of split faces
 */
CUDA_DH size_t getNumSplitFaces()
{
  return m_NumSplitFaces;
}

/**
 * @brief Get the pointer to the split faces field
 * @return the pointer to the split faces field
 */
CUDA_DH splitface_t* getSplitFaces()
{
  return m_SplitFaces;
}

CUDA_DH bool* getIsInsideCell()
{
  return m_IsInsideCell;
}

CUDA_DH bool isInsideCell(size_t idx)
{
  return m_IsInsideCell[idx];
}

CUDA_DH bool* getIsSplitCell()
{
  return m_IsInsideCell;
}

CUDA_DH bool isSplitCell(size_t idx)
{
  return m_IsInsideCell[idx];
}

CUDA_DH bool isSplitFace(size_t idx1, size_t idx2)
{
  return logicalXor(m_IsInsideCell[idx1], m_IsInsideCell[idx2]);
}

/**
 * Copy simple data attributes from another object.
 * The other object can have a different type as long as the required attributes are present.
 * param obj a constant reference to the other object
 */
template <typename T>
CUDA_HO void copyAttributes(T* obj)
{
  m_NumFields               = obj->numFields();
  m_NumVariables            = obj->numVariables();
  m_FieldSize               = obj->fieldSize();
  m_VariableSize            = obj->variableSize();
  m_NumDonorPatches         = obj->getNumDonorPatches();
  m_NumReceivingCellsConcat = obj->getNumReceivingCellsConcat();
  m_NumReceivingCellsUnique = obj->getNumReceivingCellsUnique();
  m_NumDonorWIConcat        = obj->getNumDonorWIConcat();
  m_PatchGrid               = obj->getPatchGrid();
  m_MyIndex                 = obj->getIndex();
  m_NumSplitFaces           = obj->getNumSplitFaces();
}



