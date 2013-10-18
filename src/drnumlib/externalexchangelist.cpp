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

#include "externalexchangelist.h"

#include <vector>
#include <algorithm>

#include "stringtools.h"
#include "patchgrid.h"


bool ExternalExchangeList::xyz_cell_t::operator <(const xyz_cell_t &cell) const
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

ExternalExchangeList::ExternalExchangeList(std::string name, int num_arrays, MpiCommunicator *mpi_comm, SharedMemory *shmem, Barrier *barrier)
{
  m_Name = name;
  m_MpiComm = mpi_comm;
  m_Finalised = false;
  m_SharedMem = shmem;
  m_Barrier = barrier;
  if (m_SharedMem && !m_Barrier) {
    BUG;
  }
  m_Data.resize(num_arrays);
}

void ExternalExchangeList::addCell(int grid, int index, real x, real y, real z)
{
  xyz_cell_t cell;
  cell.grid = grid;
  cell.index = index;
  cell.x = x;
  cell.y = y;
  cell.z = z;
  m_XyzCells.push_back(cell);
}

int ExternalExchangeList::grid(int i)
{
  if (m_Finalised) {
    return m_Cells[i].grid;
  } else {
    return m_XyzCells[i].grid;
  }
}

int ExternalExchangeList::index(int i)
{
  if (m_Finalised) {
    return m_Cells[i].index;
  } else {
    return m_XyzCells[i].index;
  }
}

real&  ExternalExchangeList::data(int i, int j)
{
  if (!m_Finalised) {
    BUG;
  }
  return m_Data[i][j];
}

real ExternalExchangeList::x(int i)
{
  if (m_Finalised) {
    BUG;
  }
  return m_XyzCells[i].x;
}

real ExternalExchangeList::y(int i)
{
  if (m_Finalised) {
    BUG;
  }
  return m_XyzCells[i].y;
}

real ExternalExchangeList::z(int i)
{
  if (m_Finalised) {
    BUG;
  }
  return m_XyzCells[i].z;
}

void ExternalExchangeList::mpiSend(int rank, int pos, int length)
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

void ExternalExchangeList::ipcSend(int pos, int length)
{
  if (!m_SharedMem) {
    BUG;
  }
  try {
    if (m_Finalised) {
      if (length < 0) {
        length = m_Cells.size();
      }
      std::vector<real> data(m_Cells.size());
      for (size_t i = 0; i < m_Data.size(); ++i) {
        std::string name = m_Name + "-data-" + StringTools::toString(i);

        // get existing array from shared memory (if it exists)
        int i_array = m_SharedMem->arrayIndex(name);
        if (i_array > 0) {
          if (m_SharedMem->arrayLength(i_array) != m_Cells.size()) {
            cout << m_SharedMem->arrayLength(i_array) << ", " << m_Cells.size() << endl;
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
      m_SharedMem->writeValue(m_Name + "-length", &N);
      if (N > 0) {
        m_SharedMem->writeArray(m_Name + "-cells", N, &cells[0]);
      }
    }
  } catch (IpcException e) {
    e.print();
  }
}

void ExternalExchangeList::mpiReceive(int rank, int pos, int length)
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
    m_XyzCells.resize(cells.size());
    std::copy(cells.begin(), cells.end(), m_XyzCells.begin());
  }
}

void ExternalExchangeList::ipcReceive(int pos, int length)
{
  if (!m_SharedMem) {
    BUG;
  }
  try {
    if (m_Finalised) {
      if (length < 0) {
        length = m_Cells.size();
      }
      std::vector<real> data(m_Cells.size());
      for (size_t i = 0; i < m_Data.size(); ++i) {
        std::string name = m_Name + "-data-" + StringTools::toString(i);
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
      m_SharedMem->readValue(m_Name + "-length", &N);
      std::vector<xyz_cell_t> cells(N);
      if (N > 0) {
        m_SharedMem->readArray(m_Name + "-cells", &cells[0]);
      }
      m_XyzCells.resize(cells.size());
      std::copy(cells.begin(), cells.end(), m_XyzCells.begin());
    }
  } catch (IpcException e) {
    e.print();
  }
}

void ExternalExchangeList::operator+=(const ExternalExchangeList &exchange_list)
{
  m_XyzCells.insert(m_XyzCells.end(), exchange_list.m_XyzCells.begin(), exchange_list.m_XyzCells.end());
}

void ExternalExchangeList::append(int grid, int index, real x, real y, real z)
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

void ExternalExchangeList::sort()
{
  BUG;
  if (m_Finalised) {
    BUG;
  }
  std::sort(m_XyzCells.begin(), m_XyzCells.end());

  // the follwing line is used to "squeeze" the vector
  // see http://stackoverflow.com/questions/253157/how-to-downsize-stdvector
  std::vector<xyz_cell_t>(m_XyzCells).swap(m_XyzCells);
}

void ExternalExchangeList::finalise(PatchGrid *patch_grid, int id_patch)
{
  if (patch_grid && id_patch >= 0) {
    vector<bool> cell_used(patch_grid->getPatch(id_patch)->variableSize(), false);
    real max_dist = 0;
    for (int i = 0; i < m_XyzCells.size(); ++i) {
      vec3_t xo_of = vec3_t(m_XyzCells[i].x, m_XyzCells[i].y, m_XyzCells[i].z);
      int id_cell = patch_grid->getPatch(id_patch)->findCell(xo_of);
      if (id_cell < 0) {
        ERROR("unable to find coupled cell");
      }
      vec3_t xo_dn = patch_grid->getPatch(id_patch)->xyzoCell(id_cell);
      max_dist = max(max_dist, (xo_of - xo_dn).abs());
      if (cell_used[id_cell]) {
        QString msg = "cell #" + QString::number(id_cell) + " has already been used!";
        ERROR(qPrintable(msg));
      }
      cell_used[id_cell] = true;
      m_XyzCells[i].index = id_cell;
      m_XyzCells[i].grid  = id_patch;
    }
    cout << "maximal match distance : " << max_dist << endl;
  } else {
    //sort();
  }

  m_Cells.resize(m_XyzCells.size());
  for (int i = 0; i < m_Cells.size(); ++i) {
    m_Cells[i].grid  = m_XyzCells[i].grid;
    m_Cells[i].index = m_XyzCells[i].index;
  }
  m_XyzCells.clear();

  for (int i_array = 0; i_array < m_Data.size(); ++i_array) {
    m_Data[i_array].resize(m_Cells.size());
  }

  m_Finalised = true;
}

int ExternalExchangeList::size()
{
  if (m_Finalised) {
    return m_Cells.size();
  } else {
    return m_XyzCells.size();
  }
}


