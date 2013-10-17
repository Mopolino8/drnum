
mpi_comm.barrier();
if (mpi_comm.rank() == 0) {
  barrier->wait();
}
mpi_comm.barrier();

// collect of2dn lists from MPI processes
{
  int end_cell = of2dn_list.size();
  if (mpi_comm.rank() == 0 && mpi_comm.size() > 1) {
    end_cell = of2dn_start[1];
  }
  for (int i_cell = 0; i_cell < end_cell; ++i_cell) {
    of2dn_list.data(0, i_cell) = p[of2dn_list.index(i_cell)];
    of2dn_list.data(1, i_cell) = U[of2dn_list.index(i_cell)][0];
    of2dn_list.data(2, i_cell) = U[of2dn_list.index(i_cell)][1];
    of2dn_list.data(3, i_cell) = U[of2dn_list.index(i_cell)][2];
    of2dn_list.data(4, i_cell) = thermo.T()[of2dn_list.index(i_cell)];
  }
  
  if (mpi_comm.rank() == 0) {
    for (int i_rank = 1; i_rank < mpi_comm.size(); ++i_rank) {
      of2dn_list.mpiReceive(i_rank, of2dn_start[i_rank], of2dn_size[i_rank]);
    }
  } else {
    of2dn_list.mpiSend(0);
  }
}

mpi_comm.barrier();

if (mpi_comm.rank() == 0) {
  of2dn_list.ipcSend();
  int write = 0;
  if (runTime.outputTime()) {
    write = 1;
  }
  shmem->writeValue("write", &write);
  real dt = runTime.deltaT().value();
  shmem->writeValue("dt", &dt);
  int stop = 0;
  shmem->writeValue("stop", &stop);
  barrier->wait();
  dn2of_list.ipcReceive();
}

mpi_comm.barrier();

// distribute dn2of lists to MPI processes
{  
  if (mpi_comm.rank() == 0) {
    for (int i_rank = 1; i_rank < mpi_comm.size(); ++i_rank) {
      dn2of_list.mpiSend(i_rank, dn2of_start[i_rank], dn2of_size[i_rank]);
    }
  } else {
    dn2of_list.mpiReceive(0);
  }
}




