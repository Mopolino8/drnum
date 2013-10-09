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

#ifndef EXTERNALEXCHANGELIST_H
#define EXTERNALEXCHANGELIST_H

#include <vector>
#include <algorithm>

#include "mpicommunicator.h"
#include "sharedmemory.h"
#include "barrier.h"
#include "stringtools.h"

class ExternalExchangeList
{


public: // data types


  struct cell_t
  {
    int grid, index;
  };

  struct xyz_cell_t
  {
    int grid, index;
    double x, y, z;

    bool operator<(const xyz_cell_t &cell) const
    {
      if (grid < cell.grid) {
        return true;
      } else if (grid == cell.grid  &&  index < cell.index) {
        return true;
      } else if (index == cell.index  &&  x < cell.x) {
        return true;
      } else if (x == cell.x  &&  y < cell.y) {
        return true;
      } else if (y == cell.y  &&  z < cell.z) {
        return true;
      } else {
        return false;
      }
    }

  };




private: // attributes

  std::vector<xyz_cell_t>          m_XyzCells;
  std::vector<cell_t>              m_Cells;
  MpiCommunicator                 *m_MpiComm;
  SharedMemory                    *m_SharedMem;
  Barrier                         *m_Barrier;
  bool                             m_Finalised;
  std::vector<std::vector<real> >  m_Data;


public:

  ExternalExchangeList(int num_arrays, MpiCommunicator *mpi_comm, SharedMemory *shmem = NULL, Barrier *barrier = NULL)
  {
    m_MpiComm = mpi_comm;
    m_Finalised = false;
    m_SharedMem = shmem;
    m_Barrier = barrier;
    if (m_SharedMem && !m_Barrier) {
      BUG;
    }
    m_Data.resize(num_arrays);
  }

  void addCell(int grid, int index, double x, double y, double z)
  {
    xyz_cell_t cell;
    cell.grid = grid;
    cell.index = index;
    cell.x = x;
    cell.y = y;
    cell.z = z;
    m_XyzCells.push_back(cell);
  }

  int grid(int i)
  {
    if (m_Finalised) {
      return m_Cells[i].grid;
    } else {
      return m_XyzCells[i].grid;
    }
  }

  int index(int i)
  {
    if (m_Finalised) {
      return m_Cells[i].index;
    } else {
      return m_XyzCells[i].index;
    }
  }

  real&  data(int i, int j)
  {
    if (!m_Finalised) {
      BUG;
    }
    return m_Data[i][j];
  }

  double x(int i)
  {
    if (m_Finalised) {
      BUG;
    }
    return m_XyzCells[i].x;
  }

  double y(int i)
  {
    if (m_Finalised) {
      BUG;
    }
    return m_XyzCells[i].y;
  }

  double z(int i)
  {
    if (m_Finalised) {
      BUG;
    }
    return m_XyzCells[i].z;
  }

  void mpiSend(int rank, int pos = 0, int length = -1)
  {
    if (!m_MpiComm) {
      BUG;
    }
    if (m_Finalised) {
      if (length < 0) {
        length = m_Cells.size();
      }
      for (size_t i = 0; i < m_Data.size(); ++i) {
        m_MpiComm->send(&m_Data[i][pos], length, rank, MpiCommunicator::blocking);
      }
    } else {
      if (pos != 0 || length != -1) {
        BUG;
      }
      std::vector<xyz_cell_t> cells(m_XyzCells.size());
      std::copy(m_XyzCells.begin(), m_XyzCells.end(), cells.begin());
      int N = cells.size();
      m_MpiComm->send(N, rank, MpiCommunicator::blocking);
      if (N > 0) {
        m_MpiComm->send(&cells[0], N, rank, MpiCommunicator::blocking);
      }
    }
  }

  void ipcSend(int pos = 0, int length = -1)
  {
    if (!m_SharedMem) {
      BUG;
    }
    if (m_Finalised) {
      if (length < 0) {
        length = m_Cells.size();
      }
      std::vector<real> data(m_Cells.size());
      for (size_t i = 0; i < m_Data.size(); ++i) {
        std::string name = "data-client-" + StringTools::toString(i);

        // get existing array from shared memory (if it exists)
        int i_array = m_SharedMem->arrayIndex(name);
        if (i_array) {
          if (m_SharedMem->arrayLength(i_array) != m_Cells.size()) {
            BUG;
          }
          m_SharedMem->readArray(name, &data[0]);
        }
        for (int i_data = pos; i_data < pos + length; ++i_data) {
          data[i_data] = m_Data[i][i_data];
        }
        m_SharedMem->writeArray(name, data.size(), &data[0]);
      }
    } else {
      if (pos != 0 || length != -1) {
        BUG;
      }
      std::vector<xyz_cell_t> cells(m_XyzCells.size());
      std::copy(m_XyzCells.begin(), m_XyzCells.end(), cells.begin());
      int N = cells.size();
      m_SharedMem->writeValue("length", &N);
      if (N > 0) {
        m_SharedMem->writeArray("cells", N, &cells[0]);
      }
    }
  }

  void mpiReceive(int rank, int pos = 0, int length = -1)
  {
    if (!m_MpiComm) {
      BUG;
    }
    if (m_Finalised) {
      if (length < 0) {
        length = m_Cells.size();
      }
      for (int i = 0; i < m_Data.size(); ++i) {
        m_MpiComm->receive(&m_Data[i][pos], length, rank, MpiCommunicator::blocking);
      }
    } else {
      if (pos != 0 || length != -1) {
        BUG;
      }
      int N;
      m_MpiComm->receive(N, rank, MpiCommunicator::blocking);
      std::vector<xyz_cell_t> cells(N);
      if (N > 0) {
        m_MpiComm->receive(&cells[0], N, rank, MpiCommunicator::blocking);
      }
      m_XyzCells.clear();
      m_XyzCells.insert(m_XyzCells.end(), cells.begin(), cells.end());
    }
  }

  void ipcReceive(int pos = 0, int length = -1)
  {
    if (!m_SharedMem) {
      BUG;
    }
    if (m_Finalised) {
      if (length < 0) {
        length = m_Cells.size();
      }
      std::vector<real> data(m_Cells.size());
      for (size_t i = 0; i < m_Data.size(); ++i) {
        std::string name = "data-drnum-" + StringTools::toString(i);
        m_SharedMem->readArray(name, &data[0]);
        for (int i_data = pos; i_data < pos + length; ++i_data) {
          m_Data[i][i_data] = data[i_data];
        }
      }
    } else {
      if (pos != 0 || length != -1) {
        BUG;
      }
      int N;
      std::vector<xyz_cell_t> cells(N);
      m_SharedMem->readValue("length", &N);
      if (N > 0) {
        m_SharedMem->readArray("cells", &cells[0]);
      }
      m_XyzCells.clear();
      m_XyzCells.insert(m_XyzCells.end(), cells.begin(), cells.end());
    }
  }

  void operator+=(const ExternalExchangeList &exchange_list)
  {
    m_XyzCells.insert(m_XyzCells.end(), exchange_list.m_XyzCells.begin(), exchange_list.m_XyzCells.end());
  }

  void append(int grid, int index, double x, double y, double z)
  {
    if (m_Finalised) {
      BUG;
    }
    xyz_cell_t cell;
    cell.grid = grid;
    cell.index = index;
    cell.x = x;
    cell.y = y;
    cell.z = z;
    m_XyzCells.push_back(cell);
  }

  void sort()
  {
    if (m_Finalised) {
      BUG;
    }
    std::sort(m_XyzCells.begin(), m_XyzCells.end());

    // the follwing line is used to "squeeze" the vector
    // see http://stackoverflow.com/questions/253157/how-to-downsize-stdvector
    std::vector<xyz_cell_t>(m_XyzCells).swap(m_XyzCells);
  }

  void finalise()
  {
    sort();
    m_Cells.resize(m_XyzCells.size());
    for (int i = 0; i < m_Cells.size(); ++i) {
      m_Cells[i].grid  = m_XyzCells[i].grid;
      m_Cells[i].index = m_XyzCells[i].index;
    }
    m_XyzCells.clear();
    m_Data.resize(m_Cells.size());
    m_Finalised = true;
  }

  int size()
  {
    if (m_Finalised) {
      return m_Cells.size();
    } else {
      return m_XyzCells.size();
    }
  }


};




#endif // EXTERNALEXCHANGELIST_H
