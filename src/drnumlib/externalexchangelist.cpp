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
#include <algorithm>

bool ExternalExchangeList::xyz_cell_t::operator<(const xyz_cell_t &cell)
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

bool ExternalExchangeList::cell_t::operator<(const cell_t &cell)
{
  if (grid < cell.grid) {
    return true;
  } else if (grid == cell.grid  &&  index < cell.index) {
    return true;
  } else {
    return false;
  }
}


ExternalExchangeList::ExternalExchangeList(MpiCommunicator *mpi_comm)
{
  m_MpiComm = mpi_comm;
  m_Finalised = false;
}

void ExternalExchangeList::addCell(int grid, int index, double x, double y, double z)
{
  xyz_cell_t cell;
  cell.grid = grid;
  cell.index = index;
  cell.x = x;
  cell.y = y;
  cell.z = z;
  m_XyzCells.push_back(cell);
}

void ExternalExchangeList::send()
{
  std::vector<xyz_cell_t> cells(m_XyzCells.size());
  std::copy(m_XyzCells.begin(), m_XyzCells.end(), cells.begin());
  int N = cells.size();
  m_MpiComm->send(N, 0, MpiCommunicator::blocking);
  if (N > 0) {
    m_MpiComm->send(&cells[0], N, 0, MpiCommunicator::blocking);
  }
}

void ExternalExchangeList::receive(int rank)
{
  int N;
  m_MpiComm->receive(N, rank, MpiCommunicator::blocking);
  std::vector<xyz_cell_t> cells(N);
  if (N > 0) {
    m_MpiComm->receive(&cells[0], N, rank, MpiCommunicator::blocking);
  }
  m_XyzCells.clear();
  m_XyzCells.insert(m_XyzCells.end(), cells.begin(), cells.end());
}

void ExternalExchangeList::operator +=(const ExternalExchangeList &exchange_list)
{
  m_XyzCells.insert(m_XyzCells.end(), exchange_list.m_XyzCells.begin(), exchange_list.m_XyzCells.end());
}

void ExternalExchangeList::append(int grid, int index, double x, double y, double z)
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
  if (m_Finalised) {
    BUG;
  }
  m_XyzCells.sort();
  m_Cells.sort();
}

void ExternalExchangeList::finalise()
{
  m_Cells.clear();
  for (list<xyz_cell_t>::iterator i = m_XyzCells.begin(); i != m_XyzCells.end(); ++i) {
    cell_t cell;
    cell.grid = i->grid;
    cell.index - i->index;
    m_Cells.push_back(cell);
  }
  m_XyzCells.clear();
  m_Cells.sort();
  m_Finalised = true;
}

void ExternalExchangeList::initIteration()
{
  if (m_Finalised) {
    m_Iter = m_Cells.begin();
  } else {
    m_XyzIter = m_XyzCells.begin();
  }
}

void ExternalExchangeList::nextCell()
{
  if (m_Finalised) {
    ++m_Iter;
  } else {
    ++m_XyzIter;
  }
}

bool ExternalExchangeList::endCell()
{
  if (m_Finalised) {
    return m_Iter == m_Cells.end();
  } else {
    return m_XyzIter == m_XyzCells.end();
  }
}

int ExternalExchangeList::grid()
{
  if (m_Finalised) {
    return m_Iter->grid;
  } else {
    return m_XyzIter->grid;
  }
}

int ExternalExchangeList::index()
{
  if (m_Finalised) {
    return m_Iter->index;
  } else {
    return m_XyzIter->index;
  }
}

double ExternalExchangeList::x()
{
  if (m_Finalised) {
    BUG;
  } else {
    return m_XyzIter->y;
  }
}

double ExternalExchangeList::y()
{
  if (m_Finalised) {
    BUG;
  } else {
    return m_XyzIter->y;
  }
}

double ExternalExchangeList::z()
{
  if (m_Finalised) {
    BUG;
  } else {
    return m_XyzIter->z;
  }
}

int ExternalExchangeList::size()
{
  if (m_Finalised) {
    return m_Cells.size();
  } else {
    return m_XyzCells.size();
  }
}

