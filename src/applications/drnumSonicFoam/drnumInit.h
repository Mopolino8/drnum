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

MpiCommunicator mpi_comm(argc, argv);

mpi_comm.barrier();

label id_patch = mesh.boundaryMesh().findPatchID("DrNUM"); 
const polyPatch& drnum_patch = mesh.boundaryMesh()[id_patch]; 

std::vector<std::list<int> > exchange_cells(4);
std::vector<bool> cell_marked(mesh.cells().size(), false);

Info << int(cell_marked.size()) << endl;
Info << drnum_patch.size() << endl;

forAll (drnum_patch.faceCells(), i_cell) {
  label id_cell = drnum_patch.faceCells()[i_cell];
  if (!cell_marked[id_cell]) {
    cell_marked[id_cell] = true;
    exchange_cells[0].push_back(id_cell);
  }
}

for (int i_layer = 1; i_layer < 4; ++i_layer) {
  for (std::list<int>::iterator i = exchange_cells[i_layer - 1].begin(); i != exchange_cells[i_layer - 1].end(); ++i) {
    forAll(mesh.cellCells(*i), i_cell) {
      label id_cell = mesh.cellCells(*i)[i_cell];       
      if (!cell_marked[id_cell]) {
        cell_marked[id_cell] = true;
        exchange_cells[i_layer].push_back(id_cell);
      }
    }
  }
} 

SharedMemory *shmem   = NULL;
Barrier      *barrier = NULL;

if (mpi_comm.rank() == 0) {
  try {
    shmem = new SharedMemory(1, 32*1024*1024, false);
    barrier = new Barrier(2, false);
  } catch (IpcException E) {
    E.print();
  }
}

ExternalExchangeList recv_list("recv", 5, &mpi_comm, shmem, barrier);
ExternalExchangeList send_list("send", 5, &mpi_comm, shmem, barrier);

for (int i_layer = 0; i_layer < 2; ++i_layer) {
  for (std::list<int>::iterator i = exchange_cells[i_layer].begin(); i != exchange_cells[i_layer].end(); ++i) {
    label id_cell = *i;
    double x = mesh.cellCentres()[id_cell][0];
    double y = mesh.cellCentres()[id_cell][1];
    double z = mesh.cellCentres()[id_cell][2];
    recv_list.append(mpi_comm.rank(), id_cell, x, y, z);
  }
}

for (int i_layer = 2; i_layer < 4; ++i_layer) {
  for (std::list<int>::iterator i = exchange_cells[i_layer].begin(); i != exchange_cells[i_layer].end(); ++i) {
    label id_cell = *i;
    double x = mesh.cellCentres()[id_cell][0];
    double y = mesh.cellCentres()[id_cell][1];
    double z = mesh.cellCentres()[id_cell][2];
    send_list.append(mpi_comm.rank(), id_cell, x, y, z);
  }
}

recv_list.sort();
send_list.sort();

if (mpi_comm.rank() == 0) {
  for (int i_rank = 1; i_rank < mpi_comm.size(); ++i_rank) {
    ExternalExchangeList add_list("tmp", 5, &mpi_comm, shmem, barrier);
    add_list.mpiReceive(i_rank);
    recv_list += add_list;
    add_list.mpiReceive(i_rank);
    send_list += add_list;
  }
} else {
  recv_list.mpiSend(0);
  send_list.mpiSend(0);
}

// some infos -- might be deleted later on

Info << endl;
Info << "receive cells:\n";
Info << "--------------" << endl;

{
  int i = 0;
  for (int i_rank = 0; i_rank < mpi_comm.size(); ++i_rank) {
    int N = 0;
    while (recv_list.grid(i) == i_rank && i < recv_list.size()) {
      ++N;
      ++i;
    }
    Info << "rank = " << i_rank << "   " << N << " cells" << endl;
  }
}

Info << endl;
Info << "send cells:\n";
Info << "-----------" << endl;

{
  int i = 0;
  for (int i_rank = 0; i_rank < mpi_comm.size(); ++i_rank) {
    int N = 0;
    while (send_list.grid(i) == i_rank && i < send_list.size()) {
      ++N;
      ++i;
    }
    Info << "rank = " << i_rank << "   " << N << " cells" << endl;
  }
}

// end of infos

if (mpi_comm.rank() == 0) {
  Info << "trying to connect to DrNUM ..." << endl;
  int client_ready = 1;
  shmem->writeValue("client-ready", &client_ready);
  barrier->wait();
  Info << "connection established" << endl;

  send_list.ipcSend();
  recv_list.ipcSend();
  barrier->wait();
  barrier->wait();
}

Info << "exiting" << endl;
exit(-1);



