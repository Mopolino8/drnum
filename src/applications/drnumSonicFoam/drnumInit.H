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

ExternalExchangeList dn2of_list("dn2of", 5, &mpi_comm, shmem, barrier);
ExternalExchangeList of2dn_list("of2dn", 5, &mpi_comm, shmem, barrier);

for (int i_layer = 0; i_layer < 2; ++i_layer) {
  for (std::list<int>::iterator i = exchange_cells[i_layer].begin(); i != exchange_cells[i_layer].end(); ++i) {
    label id_cell = *i;
    double x = mesh.cellCentres()[id_cell][0];
    double y = mesh.cellCentres()[id_cell][1];
    double z = mesh.cellCentres()[id_cell][2];
    dn2of_list.append(mpi_comm.rank(), id_cell, x, y, z);
  }
}

for (int i_layer = 2; i_layer < 4; ++i_layer) {
  for (std::list<int>::iterator i = exchange_cells[i_layer].begin(); i != exchange_cells[i_layer].end(); ++i) {
    label id_cell = *i;
    double x = mesh.cellCentres()[id_cell][0];
    double y = mesh.cellCentres()[id_cell][1];
    double z = mesh.cellCentres()[id_cell][2];
    of2dn_list.append(mpi_comm.rank(), id_cell, x, y, z);
  }
}

dn2of_list.sort();
of2dn_list.sort();

// create label lists in order to fix the cells during computation
labelList fixed_cells(dn2of_list.size());

if (mpi_comm.rank() == 0) {
  for (int i_rank = 1; i_rank < mpi_comm.size(); ++i_rank) {
    ExternalExchangeList add_list("tmp", 5, &mpi_comm, shmem, barrier);
    add_list.mpiReceive(i_rank);
    dn2of_list += add_list;
    add_list.mpiReceive(i_rank);
    of2dn_list += add_list;
  }
} else {
  dn2of_list.mpiSend(0);
  of2dn_list.mpiSend(0);
}

std::vector<int> of2dn_start(mpi_comm.size(), 0);
std::vector<int> dn2of_start(mpi_comm.size(), 0);
std::vector<int> of2dn_size(mpi_comm.size(), 0);
std::vector<int> dn2of_size(mpi_comm.size(), 0);

if (mpi_comm.rank() == 0) {
  
  Info << endl;
  Info << "dn2of cells:\n";
  Info << "----- ------" << endl;

  {
    int i = 0;
    int N_total = 0;
    for (int i_rank = 0; i_rank < mpi_comm.size(); ++i_rank) {
      int N = 0;
      dn2of_start[i_rank] = i;
      while (dn2of_list.grid(i) == i_rank && i < dn2of_list.size()) {
        ++N;
        ++dn2of_size[i_rank];
        ++i;
      }
      Info << "rank  = " << i_rank << "   " << N << " cells" << endl;
      N_total += N;
    }
    Info << "total = " << N_total << " cells" << endl;
  }

  Info << endl;
  Info << "of2dn cells:\n";
  Info << "------------" << endl;

  {
    int i = 0;
    int N_total = 0;
    for (int i_rank = 0; i_rank < mpi_comm.size(); ++i_rank) {
      int N = 0;
      of2dn_start[i_rank] = i;
      while (of2dn_list.grid(i) == i_rank && i < of2dn_list.size()) {
        ++N;
        ++of2dn_size[i_rank];
        ++i;
      }
      Info << "rank  = " << i_rank << "   " << N << " cells" << endl;
      N_total += N;
    }
    Info << "total = " << N_total << " cells" << endl;
  }

  Info << "\ntrying to connect to DrNUM ..." << endl;
  int client_ready = 1;
  shmem->writeValue("client-ready", &client_ready);
  barrier->wait();
  Info << "connection established" << endl;

  of2dn_list.ipcSend();
  dn2of_list.ipcSend();
  barrier->wait();
  barrier->wait();
}

of2dn_list.finalise();
dn2of_list.finalise();




