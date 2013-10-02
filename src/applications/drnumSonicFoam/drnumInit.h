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

forAll(drnum_patch, id_face) {
  label id_cell = mesh.owner()[id_face];
  if (!cell_marked[id_cell]) {
    cell_marked[id_cell] = true;
    exchange_cells[0].push_back(id_cell);
  }
}

Info << int(exchange_cells[0].size()) << endl;
for (int i_layer = 1; i_layer < 4; ++i_layer) {
  for (std::list<int>::iterator i = exchange_cells[i_layer - 1].begin(); i != exchange_cells[i_layer - 1].end(); ++i) {
    //Info << i_layer << "," << *i << ": " << mesh.cellCells(*i).size() << endl;
    forAll(mesh.cellCells(*i), id_cell) {
      if (!cell_marked[id_cell]) {
        cell_marked[id_cell] = true;
        exchange_cells[i_layer].push_back(id_cell);
      }
    }
  }
  Info << int(exchange_cells[i_layer].size()) << endl;
} 

ExternalExchangeList receive_list(&mpi_comm);
ExternalExchangeList send_list(&mpi_comm);

for (int i_layer = 0; i_layer < 2; ++i_layer) {
  for (std::list<int>::iterator i = exchange_cells[i_layer].begin(); i != exchange_cells[i_layer].end(); ++i) {
    label id_cell = *i;
    double x = mesh.cellCentres()[id_cell][0];
    double y = mesh.cellCentres()[id_cell][1];
    double z = mesh.cellCentres()[id_cell][2];
    receive_list.append(mpi_comm.rank(), id_cell, x, y, z);
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

receive_list.sort();
send_list.sort();

if (mpi_comm.rank() == 0) {
  for (int i_rank = 1; i_rank < mpi_comm.size(); ++i_rank) {
    ExternalExchangeList add_list(&mpi_comm);
    add_list.receive(i_rank);
    receive_list += add_list;
    add_list.receive(i_rank);
    send_list += add_list;
  }
} else {
  receive_list.send();
  send_list.send();
}

Info << endl;
Info << "receive cells:\n";
Info << "--------------" << endl;

receive_list.initIteration();
for (int i_rank = 0; i_rank < mpi_comm.size(); ++i_rank) {
  int N = 0;
  while (receive_list.grid() == i_rank && !receive_list.endCell()) {
    ++N;
    receive_list.nextCell();
  }
  Info << "rank = " << i_rank << "   " << N << " cells" << endl;
}

Info << endl;
Info << "send cells:\n";
Info << "-----------" << endl;

send_list.initIteration();
for (int i_rank = 0; i_rank < mpi_comm.size(); ++i_rank) {
  int N = 0;
  while (send_list.grid() == i_rank && !send_list.endCell()) {
    ++N;
    send_list.nextCell();
  }
  Info << "rank = " << i_rank << "   " << N << " cells" << endl;
}

mpi_comm.barrier();
Info << "exiting" << endl;
exit(-1);




//Info << "Area " << A << endl; 
