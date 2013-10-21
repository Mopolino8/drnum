{
  int start = 0;
  int stop = dn2of_list.size();
  if (mpi_comm.rank() == 0) {
    start = dn2of_start[0];
    stop = dn2of_start[0] + dn2of_size[0];
  }
  for (int i_cell = start; i_cell < stop; ++i_cell) {
    label cell = dn2of_list.index(i_cell);
    scalar p_dn = dn2of_list.data(0, i_cell);
    scalar u_dn = dn2of_list.data(1, i_cell);
    scalar v_dn = dn2of_list.data(2, i_cell);
    scalar w_dn = dn2of_list.data(3, i_cell);
    scalar T_dn = dn2of_list.data(4, i_cell);
    U[cell][0] = u_dn;
    U[cell][1] = v_dn;
    U[cell][2] = w_dn;
    scalar e_of = e[cell];
    scalar T_of = thermo.T()[cell];
    scalar cv   = e_of/T_of;
    scalar e_dn = cv*T_dn;
    p[cell] = p_dn;
    e[cell] = e_dn;
  }
  thermo.correct();
}



