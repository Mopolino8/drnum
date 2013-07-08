#ifndef STRUCTUREDHEXRASTER_H
#define STRUCTUREDHEXRASTER_H

// #include <cstddef>
#include "raster.h"

//class StructuredHexRaster;

//#ifdef WITH_VTK
//#include <QString>
//#include <QVector>
//#include <vtkSmartPointer.h>
//#include <vtkRectilinearGrid.h>
//#include <vtkXMLRectilinearGridWriter.h>
//#include <vtkFloatArray.h>
//#include <vtkCellData.h>
//#endif

class StructuredHexRaster : public Raster
{

protected: // attributes

  size_t m_NumI;  ///< Number of cells aligned in i-direction
  size_t m_NumJ;  ///< Number of cells aligned in j-direction
  size_t m_NumK;  ///< Number of cells aligned in k-direction

protected: // methods

public: // methods

  StructuredHexRaster();

  virtual void resize(const size_t& num_i, const size_t& num_j, const size_t& num_k);

  size_t sizeI() { return m_NumI; }
  size_t sizeJ() { return m_NumJ; }
  size_t sizeK() { return m_NumK; }

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @return the index in the one dimensional data field
   */

  /// @todo check, why "int" and not "size_t"
  size_t index(const int& i, const int& j, const int& k) const
  {
    return i*m_NumJ*m_NumK + j*m_NumK + k;
  }

  /**
   * @brief Get the indicees (i, j, k) from field index/
   * @param l the field index
   * @param i first  index  (return reference)
   * @param j second  index (return reference)
   * @param k third  index  (return reference)
   * @return the index in the one dimensional data field
   */
  inline void ijk(const size_t& l,
                  size_t& i, size_t& j, size_t& k) const;

  /**
   * @brief Get the field index of an (i, j, k) triple/
   * @param i first  index
   * @param j second  index
   * @param k third  index
   * @param error true, if (i, j, k) are out of bounds
   * @return the index in the one dimensional data field
   */
  size_t save_index(const int& i, const int& j, const int& k,
                    bool& error) const ;

  /**
   * Check if an (i,j,k) triple is inside the patch.
   * Attention only the upper limit will be checked (unsigned data type).
   * @param i first index
   * @param j second index
   * @param k third index
   * @return true if it is a valid (i,j,k) triple, false otherwise
   */
  bool checkRange(const size_t& i, const size_t& j, const size_t& k);

};

inline bool StructuredHexRaster::checkRange(const size_t& i, const size_t& j, const size_t& k)
{
  if (i >= sizeI() || j >= sizeJ() || k >= sizeK()) {
    return false;
  }
  return true;
}


inline void StructuredHexRaster::ijk(const size_t& l,
                                     size_t& i, size_t& j, size_t& k) const
{
  size_t rest;
  i = l / (m_NumJ*m_NumK);
  rest = l - i*(m_NumJ*m_NumK);
  j = rest / m_NumK;
  k = rest - j*m_NumK;
}


inline size_t StructuredHexRaster::save_index(const int& i, const int& j, const int& k,
                                              bool& error) const
{
  size_t si = i;  // avoid vast compiler warnings
  size_t sj = j;
  size_t sk = k;
  error = false;
  if(i < 0) error = true;
  if(j < 0) error = true;
  if(k < 0) error = true;
  if(si > m_NumI-1) error = true;
  if(sj > m_NumJ-1) error = true;
  if(sk > m_NumK-1) error = true;
  return si*m_NumJ*m_NumK + sj*m_NumK + sk;
}


#endif // STRUCTUREDHEXRASTER_H
