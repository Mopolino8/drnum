#include "cubeincartisianpatch.h"
#include "perfectgas.h"

CubeInCartisianPatch::CubeInCartisianPatch()
{
  m_Patch = NULL;
  m_Start.i = 0;
  m_Start.j = 0;
  m_Start.k = 0;
  m_Stop.i = 0;
  m_Stop.j = 0;
  m_Stop.k = 0;
}

CubeInCartisianPatch::CubeInCartisianPatch(CartesianPatch *patch)
{
  m_Patch = patch;
  m_Start.i = patch->sizeI()/4;
  m_Start.j = patch->sizeJ()/4;
  m_Start.k = patch->sizeK()/4;
  m_Stop.i = 3*patch->sizeI()/4;
  m_Stop.j = 3*patch->sizeJ()/4;
  m_Stop.k = 3*patch->sizeK()/4;
  m_NumLayers = 2;
  m_Count.resize(m_Patch->fieldSize());
}



void CubeInCartisianPatch::setRange(size3_t i_start, size3_t i_stop)
{
  m_Start = i_start;
  m_Stop = i_stop;
}

void CubeInCartisianPatch::setRange(vec3_t x1, vec3_t x2)
{
  m_Start.i = m_Patch->sizeI();
  m_Start.j = m_Patch->sizeJ();
  m_Start.k = m_Patch->sizeK();
  m_Stop.i = 0;
  m_Stop.j = 0;
  m_Stop.k = 0;
  for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
        vec3_t x;
        x[0] = (real(i) + 0.5)*m_Patch->dx();
        x[1] = (real(j) + 0.5)*m_Patch->dy();
        x[2] = (real(k) + 0.5)*m_Patch->dz();
        x = m_Patch->getTransformInertial2This().transformReverse(x);
        if (x[0] > x1[0] && x[1] > x1[1] && x[2] > x1[2]) {
          if (x[0] < x2[0] && x[1] < x2[1] && x[2] < x2[2]) {
            m_Start.i = min(m_Start.i, i);
            m_Start.j = min(m_Start.j, j);
            m_Start.k = min(m_Start.k, k);
            m_Stop.i = max(m_Stop.i, i);
            m_Stop.j = max(m_Stop.j, j);
            m_Stop.k = max(m_Stop.k, k);
          }
        }
      }
    }
  }
}

void CubeInCartisianPatch::operator ()()
{
  m_Patch->copyFieldToHost(0);

  real p0 = 0;
  real T0 = 0;
  size_t N = 0;
  real var1[5], var2[5], p, T, u, v, w;

  for (int i_var = 0; i_var < 5; ++i_var) {
    var1[i_var] = 0;
  }
  for (size_t i = m_Start.i; i < m_Stop.i; ++i) {
    for (size_t j = m_Start.j; j < m_Stop.j; ++j) {
      for (size_t k = m_Start.k; k < m_Stop.k; ++k) {
        m_Patch->setVar(0, i, j, k, var1);
      }
    }
  }

  for (size_t layer = 0; layer < m_NumLayers; ++layer) {

    resetCount();

    size3_t start = m_Start;
    size3_t stop = m_Stop;
    if (m_Start.i > 0) {
      start.i += layer;
    }
    if (m_Start.j > 0) {
      start.j += layer;
    }
    if (m_Start.k > 0) {
      start.k += layer;
    }
    if (m_Stop.i < m_Patch->sizeI() - 1) {
      stop.i -= layer;
    }
    if (m_Stop.j < m_Patch->sizeJ() - 1) {
      stop.j -= layer;
    }
    if (m_Stop.k < m_Patch->sizeK() - 1) {
      stop.k -= layer;
    }

    // K planes
    for (size_t i = start.i; i < stop.i; ++i) {
      for (size_t j = start.j; j < stop.j; ++j) {
        if (m_Start.k > 0) {
          m_Patch->getVar(0, i, j, start.k - 1, var1);
          m_Patch->getVar(0, i, j, start.k, var2);
          PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);
          PerfectGas::primitiveToConservative(p, T, 0, 0, 0, var1);
          for (int i_var = 0; i_var < 5; ++i_var) {
            var2[i_var] += var1[i_var];
          }
          count(i, j, start.k);
          m_Patch->setVar(0, i, j, start.k, var2);
          if (layer == 0) {
            ++N;
            p0 += p;
            T0 += T;
          }
        }

        if (m_Stop.k < m_Patch->sizeK() - 1) {
          m_Patch->getVar(0, i, j, stop.k, var1);
          m_Patch->getVar(0, i, j, stop.k - 1, var2);
          PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);
          PerfectGas::primitiveToConservative(p, T, 0, 0, 0, var1);
          for (int i_var = 0; i_var < 5; ++i_var) {
            var2[i_var] += var1[i_var];
          }
          count(i, j, stop.k - 1);
          m_Patch->setVar(0, i, j, stop.k - 1, var2);
          if (layer == 0) {
            ++N;
            p0 += p;
            T0 += T;
          }
        }
      }
    }

    // I planes
    for (size_t j = start.j; j < stop.j; ++j) {
      for (size_t k = start.k; k < stop.k; ++k) {
        if (m_Start.i > 0) {
          m_Patch->getVar(0, start.i - 1, j, k, var1);
          m_Patch->getVar(0, start.i, j, k, var2);
          PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);
          PerfectGas::primitiveToConservative(p, T, 0, 0, 0, var1);
          for (int i_var = 0; i_var < 5; ++i_var) {
            var2[i_var] += var1[i_var];
          }
          count(start.i, j, k);
          m_Patch->setVar(0, start.i, j, k, var2);
          if (layer == 0) {
            ++N;
            p0 += p;
            T0 += T;
          }
        }

        if (m_Stop.i < m_Patch->sizeI() - 1) {
          m_Patch->getVar(0, stop.i, j, k, var1);
          m_Patch->getVar(0, stop.i - 1, j, k, var2);
          PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);
          PerfectGas::primitiveToConservative(p, T, 0, 0, 0, var1);
          for (int i_var = 0; i_var < 5; ++i_var) {
            var2[i_var] += var1[i_var];
          }
          count(stop.i - 1, j, k);
          m_Patch->setVar(0, stop.i - 1, j, k, var2);
          if (layer == 0) {
            ++N;
            p0 += p;
            T0 += T;
          }
        }
      }
    }

    // J planes
    for (size_t i = start.i; i < stop.i; ++i) {
      for (size_t k = start.k; k < stop.k; ++k) {
        if (m_Start.j > 0) {
          m_Patch->getVar(0, i, start.j - 1, k, var1);
          m_Patch->getVar(0, i, start.j, k, var2);
          PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);
          PerfectGas::primitiveToConservative(p, T, 0, 0, 0, var1);
          for (int i_var = 0; i_var < 5; ++i_var) {
            var2[i_var] += var1[i_var];
          }
          count(i, start.j, k);
          m_Patch->setVar(0, i, start.j, k, var2);
          if (layer == 0) {
            ++N;
            p0 += p;
            T0 += T;
          }
        }

        if (m_Stop.j < m_Patch->sizeJ() - 1) {
          m_Patch->getVar(0, i, stop.j, k, var1);
          m_Patch->getVar(0, i, stop.j - 1, k, var2);
          PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);
          PerfectGas::primitiveToConservative(p, T, 0, 0, 0, var1);
          for (int i_var = 0; i_var < 5; ++i_var) {
            var2[i_var] += var1[i_var];
          }
          count(i, stop.j - 1, k);
          m_Patch->setVar(0, i, stop.j - 1, k, var2);
          if (layer == 0) {
            ++N;
            p0 += p;
            T0 += T;
          }
        }
      }
    }

    for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
      for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
        for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
          if (getCount(i, j, k) > 0) {
            m_Patch->getVar(0, i, j, k, var1);
            for (int i_var = 0; i_var < 5; ++i_var) {
              var1[i_var] /= getCount(i, j, k);
            }
            m_Patch->setVar(0, i, j, k, var1);
          }
        }
      }
    }

  }

  p0 /= N;
  T0 /= N;
  {
    size3_t start = m_Start;
    size3_t stop = m_Stop;
    if (m_Start.i > 0) {
      start.i += m_NumLayers;
    }
    if (m_Start.j > 0) {
      start.j += m_NumLayers;
    }
    if (m_Start.k > 0) {
      start.k += m_NumLayers;
    }
    if (m_Stop.i < m_Patch->sizeI() - 1) {
      stop.i -= m_NumLayers;
    }
    if (m_Stop.j < m_Patch->sizeJ() - 1) {
      stop.j -= m_NumLayers;
    }
    if (m_Stop.k < m_Patch->sizeK() - 1) {
      stop.k -= m_NumLayers;
    }
    PerfectGas::primitiveToConservative(p0, T0, 0, 0, 0, var1);
    for (size_t i = start.i; i < stop.i; ++i) {
      for (size_t j = start.j; j < stop.j; ++j) {
        for (size_t k = start.k; k < stop.k; ++k) {
          m_Patch->setVar(0, i, j, k, var1);
        }
      }
    }
  }

  m_Patch->copyFieldToDevice(0);
}
