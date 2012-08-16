#ifndef TINSECTHASHRASTER_HH
#define TINSECTHASHRASTER_HH

#include <cstddef>
#include <namespace_mouse.hh>

BEGIN_MOUSE
template <class T> class TInsectHashRaster;
END_MOUSE

#include <TInsectHashRaster.hh>
#include "CartesianRaster.hh"
//#include <vector>

BEGIN_MOUSE

/** Class to hold a TInsectHashRaster on a CartesianRaster
 *
 *  Hold arbitrary ammounts of data types T on each cell of the CartesianRaster
 *
 */
template <class T>
class TInsectHashRaster : public CartesianRaster
{
  protected:
  /**
    * TInsectionList to hold data.
    */
  TInsectionList<T>* m_tinsect;


  public:
  /** constructor
    */
  TInsectHashRaster(real m_x_min, real m_y_min, real m_z_min,
                    real m_x_max, real m_y_max, real m_z_max,
                    size_t num_i, size_t num_j, size_t num_k,
                    size_t content_size = 0, size_t content_increment = 0);

  /** Initializer for contents of m_tinsect.
    *  Note: InitInsect() is automatically done in constructor. Use for reinitialization.
    */
  void initInsect();


  /**
    * Insert an element at position (i,j,k)
    * @param i i-index
    * @param j j-index
    * @param k k-index
    * @param value value to insert
    */
  void insert(size_t i, size_t j, size_t k, T value);


  /** See basesrc/templates/TInsectionList.hh
    *  Get the number of entries in the second dimension of t_insect.
    *  @param l the index of the cell in the raster
    *  @return the number of items stored for l
    */
  size_t getNumItems(size_t l);


  /** See TInsectionList.hh
    *  Data access method (read/write)
    *  @param l the index of the node in the raster
    *  @param l_i the l_i-th entry for a raster node "l"
    *  @return a reference to the entry
    */
  T& at(size_t& l, size_t& l_i);


  /** Data access method (read only)
    *  Hand over to t_insect contained.
    *  @param m the index of the node in the raster
    *  @param item_start_ptr pointer to the location in item_list, where entries of "i"
    *         start (reference return value).
    *  @param num_items number of items for entries of "i" (reference return value).
    */
  void getItemPtrAndNum(size_t l,
                        T*& item_start_ptr, size_t& num_items);


  /** destructor */
  virtual ~TInsectHashRaster() {}  /// @todo destructor missing

};

template <class T>
inline void TInsectHashRaster<T>::InitInsect()
{
  m_tinsect->Initialize();
}

template <class T>
inline void TInsectHashRaster<T>::insert(size_t i, size_t j, size_t k, T value)
{
  m_tinsect->ConsiderStore(index(i, j, k), value);
}

template <class T>
inline void TInsectHashRaster<T>::Finalize()
{
  m_tinsect->Finalize();
}

template <class T>
inline bool TInsectHashRaster<T>::save_insert(size_t i, size_t j, size_t k, T value)
{
  bool error = false;
  size_t l_index = save_index(i, j, k,
                              error);
  if (!error) {
    m_tinsect->ConsiderStore(l_index, value);
  }
  return !error;
}

template <class T>
inline size_t TInsectHashRaster<T>::getNumItems(size_t l) {
  return m_tinsect->NumItems(l);
}

template <class T>
inline T& TInsectHashRaster<T>::at(size_t& l, size_t& l_i) {
  return m_tinsect->At(l, l_i);
}

template <class T>
inline void TInsectHashRaster<T>::getItemPtrAndNum(size_t l,
                                                   T*& item_start_ptr, size_t& num_items) {
  m_tinsect->ItemPtrAndNum(l,
                           item_start_ptr, num_items);
}

#include "TInsectHashRaster.cc"
END_MOUSE

#endif  //TINSECTHASHRASTER_HH
