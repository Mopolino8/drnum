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

#include "mpicommunicator.h"
#include "sharedmemory.h"
#include "barrier.h"

class PatchGrid;

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
    real x, y, z;

    bool operator<(const xyz_cell_t &cell) const;
  };


private: // attributes

  std::vector<xyz_cell_t>          m_XyzCells;
  std::vector<cell_t>              m_Cells;
  MpiCommunicator                 *m_MpiComm;
  SharedMemory                    *m_SharedMem;
  Barrier                         *m_Barrier;
  bool                             m_Finalised;
  std::vector<std::vector<real> >  m_Data;
  std::string                      m_Name;


public:

  ExternalExchangeList(std::string name, int num_arrays, MpiCommunicator *mpi_comm, SharedMemory *shmem = NULL, Barrier *barrier = NULL);
  void addCell(int grid, int index, real x, real y, real z);
  int grid(int i);
  int index(int i);
  real& data(int i, int j);
  real x(int i);
  real y(int i);
  real z(int i);
  void mpiSend(int rank, int pos = 0, int length = -1);
  void ipcSend(int pos = 0, int length = -1);
  void mpiReceive(int rank, int pos = 0, int length = -1);
  void ipcReceive(int pos = 0, int length = -1);
  void operator+=(const ExternalExchangeList &exchange_list);
  void append(int grid, int index, real x, real y, real z);
  void sort();
  void finalise(PatchGrid *patch_grid = NULL, int id_patch = -1);
  int size();

};




#endif // EXTERNALEXCHANGELIST_H
