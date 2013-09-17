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

#ifndef MPICOMMUNICATOR_H
#define MPICOMMUNICATOR_H

#include <mpi.h>
#include <vector>

class MpiCommunicator
{

private: // attributes

  int m_Size;
  int m_Rank;

  std::vector<MPI_Request> m_SendReq;
  std::vector<bool>        m_SendPending;
  std::vector<MPI_Request> m_RecvReq;
  std::vector<bool>        m_RecvPending;


public:

  MpiCommunicator(int argc, char **argv);
  ~MpiCommunicator();

  int  size()    { return m_Size; }
  int  rank()    { return m_Rank; }
  void barrier() { MPI_Barrier(MPI_COMM_WORLD); }

  template <typename T> void broadcast(T *t, int bcast_rank, int n);
  template <typename T> void broadcast(T &t, int bcast_rank);
  template <typename T> void send(T &t, int to_rank, bool non_block);
  template <typename T> void send(T *t, int length, int to_rank, bool non_block);
  template <typename T> T    receive(int from_rank, bool non_block);
  template <typename T> void receive(T& t, int from_rank, bool non_block);
  template <typename T> void receive(T* t, int length, int from_rank, bool non_block);

};

template <typename T>
void MpiCommunicator::broadcast(T *t, int bcast_rank, int n)
{
  MPI_Bcast(t, n*sizeof(T), MPI_CHAR, bcast_rank, MPI_COMM_WORLD);
}

template <typename T>
void MpiCommunicator::broadcast(T &t, int bcast_rank)
{
  MPI_Bcast(&t, sizeof(T), MPI_CHAR, bcast_rank, MPI_COMM_WORLD);
}

template <typename T>
void MpiCommunicator::send(T &t, int to_rank, bool non_block)
{
  MPI_Status status;
  if (m_SendPending[to_rank]) MPI_Wait(&m_SendReq[to_rank], &status);
  int tag = 1000*rank() + to_rank;
  MPI_Isend(&t, sizeof(T), MPI_CHAR, to_rank, tag, MPI_COMM_WORLD, &m_SendReq[to_rank]);
  if (non_block) {
    m_SendPending[to_rank] = true;
  } else {
    MPI_Wait(&m_SendReq[to_rank], &status);
    m_SendPending[to_rank] = false;
  }
}

template <typename T>
void MpiCommunicator::send(T *t, int length, int to_rank, bool non_block)
{
  MPI_Status status;
  if (m_SendPending[to_rank]) MPI_Wait(&m_SendReq[to_rank], &status);
  int tag = 1000*rank() + to_rank;
  MPI_Isend(t, sizeof(T)*length, MPI_CHAR, to_rank, tag, MPI_COMM_WORLD, &m_SendReq[to_rank]);
  if (non_block) {
    m_SendPending[to_rank] = true;
  } else {
    MPI_Wait(&m_SendReq[to_rank], &status);
    m_SendPending[to_rank] = false;
  }
}

template <typename T>
T MpiCommunicator::receive(int from_rank, bool non_block)
{
  MPI_Status status;
  if (m_RecvPending[from_rank]) MPI_Wait(&m_RecvReq[from_rank], &status);
  int tag = 1000*from_rank + rank();
  T t;
  MPI_Irecv(&t, sizeof(T), MPI_CHAR, from_rank, tag, MPI_COMM_WORLD, &m_RecvReq[from_rank]);
  if (non_block) {
    m_RecvPending[from_rank] = true;
  } else {
    MPI_Wait(&m_RecvReq[from_rank], &status);
    m_RecvPending[from_rank] = false;
  }
  return t;
}

template <typename T>
void MpiCommunicator::receive(T& t, int from_rank, bool non_block)
{
  t = receive<T>(from_rank);
}

template <typename T>
void MpiCommunicator::receive(T* t, int length, int from_rank, bool non_block)
{
  MPI_Status status;
  if (m_RecvPending[from_rank]) MPI_Wait(&m_RecvReq[from_rank], &status);
  int tag = 1000*from_rank + rank();
  MPI_Irecv(t, sizeof(T)*length, MPI_CHAR, from_rank, tag, MPI_COMM_WORLD, &m_RecvReq[from_rank]);
  if (non_block) {
    m_RecvPending[from_rank] = true;
  } else {
    MPI_Wait(&m_RecvReq[from_rank], &status);
    m_RecvPending[from_rank] = false;
  }
}


#endif // MPICOMMUNICATOR_H
