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

#include "mpicommunicator.h"

MpiCommunicator::MpiCommunicator(int argc, char **argv)
{
  int initialised;
  MPI_Initialized(&initialised);
  if (!initialised) {
    MPI_Init(&argc, &argv);
  }
  MPI_Comm_size(MPI_COMM_WORLD, &m_Size);
  MPI_Comm_rank(MPI_COMM_WORLD, &m_Rank);
  m_SendReq.resize(size());
  m_RecvReq.resize(size());
  m_SendPending.resize(size(),false);
  m_RecvPending.resize(size(),false);
}

MpiCommunicator::~MpiCommunicator()
{
  barrier();
  MPI_Finalize();
}
