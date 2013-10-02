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
#include <list>

#include "mpicommunicator.h""

class ExternalExchangeList
{


public: // data types

  struct cell_t
  {
    int grid, index;
    bool operator<(const cell_t &cell);
  };

  struct xyz_cell_t
  {
    int grid, index;
    double x, y, z;
    bool operator<(const xyz_cell_t &cell);
  };


private: // attributes

  std::list<xyz_cell_t> m_XyzCells;
  std::list<cell_t>     m_Cells;
  MpiCommunicator      *m_MpiComm;
  bool                  m_Finalised;

  std::list<xyz_cell_t>::iterator m_XyzIter;
  std::list<cell_t>::iterator     m_Iter;


public:

  ExternalExchangeList(MpiCommunicator *mpi_comm);

  void   addCell(int grid, int index, double x, double y, double z);
  int    grid();
  int    index();
  double x();
  double y();
  double z();
  void   send();
  void   receive(int rank);
  void   operator+=(const ExternalExchangeList &exchange_list);
  void   append(int grid, int index, double x, double y, double z);
  void   sort();
  void   finalise();
  int    size();
  void   initIteration();
  void   nextCell();
  bool   endCell();

};

#endif // EXTERNALEXCHANGELIST_H
