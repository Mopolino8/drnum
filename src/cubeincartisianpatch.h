#ifndef CUBEINCARTISIANPATCH_H
#define CUBEINCARTISIANPATCH_H

#include "genericoperation.h"
#include "cartesianpatch.h"

class CubeInCartisianPatch : public GenericOperation
{

  CartesianPatch *m_Patch;
  size3_t         m_Start;
  size3_t         m_Stop;
  size_t          m_NumLayers;
  vector<int>     m_Count;


  void resetCount();
  void count(size_t i, size_t j, size_t k);
  int  getCount(size_t i, size_t j, size_t k);

public:

  CubeInCartisianPatch();
  CubeInCartisianPatch(CartesianPatch* patch);
  void setRange(size3_t i_start, size3_t i_stop);
  void setNumLayers(size_t N) { m_NumLayers = N; }

  virtual void operator()();

};


inline void CubeInCartisianPatch::resetCount()
{
  for (size_t i = 0; i < m_Patch->fieldSize(); ++i) {
    m_Count[i] = 0;
  }
}

inline void CubeInCartisianPatch::count(size_t i, size_t j, size_t k)
{
  size_t idx = m_Patch->index(i, j, k);
  ++m_Count[idx];
}

inline int CubeInCartisianPatch::getCount(size_t i, size_t j, size_t k)
{
  size_t idx = m_Patch->index(i, j, k);
  return m_Count[idx];
}

#endif // CUBEINCARTISIANPATCH_H
