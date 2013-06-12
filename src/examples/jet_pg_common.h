// Original. Taken over from jet_common.h

// size_t N = 20;  => see grid files
real Ma_jet         = 0.30;
real Ma_far         = 0.15;
real p              = 1.0e5;
real p_jet          = 1.0e5;
real T              = 300;
real T_jet          = 800;
real u_jet          = Ma_jet*sqrt(PerfectGas::gamma()*PerfectGas::R()*T_jet);
real u_far          = Ma_far*sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
real Re             = 1e5;
real L              = Re*PerfectGas::mu()*PerfectGas::R()*T_jet/(p_jet*u_jet);
// real ar             = 3.0;
// real length_factor  = 2.0;
real time           = L/u_jet;
real cfl_target     = 0.2;
real t_write        = 0;
real write_interval = 0.5*time;
//real write_interval = 0.2*time;
real total_time     = 100*time;

cout << "L = " << L << endl;

// Modified bejond this line to suit PatchGrids
// Patch grid
PatchGrid patch_grid;
//.. general settings (apply to all subsequent patches)
patch_grid.setNumberOfFields(3);
patch_grid.setNumberOfVariables(5);
patch_grid.setInterpolateData();
patch_grid.setNumOverlapLayers(2);
patch_grid.setTransferPadded(true);
//.. read patch grid from file
patch_grid.readGrid(grid_file_name);
//.. re-scale grid to suit Re number needs

patch_grid.scaleRefParental(L);


//.. compute dependencies
patch_grid.computeDependencies(true);

//CartesianPatch patch;
//patch.setNumberOfFields(3);
//patch.setNumberOfVariables(5);
//patch.setupAligned(0, 0, 0, length_factor*18*L, 3*L, 3*L);
//size_t NI = 2*length_factor*3.0/ar*N;
//size_t NJ = N;
//size_t NK = N;
//patch.resize(NI, NJ, NK);

// Time step
real ch_speed = (max(u_jet, u_far) + sqrt(PerfectGas::gamma()*PerfectGas::R()*max(T_jet,T)));
//real dt = cfl_target*patch_grid.baseTimeStep(ch_speed);

real dt = cfl_target*patch_grid.computeMinChLength()/ch_speed;

cout << " patch_grid.computeMinChLength() = " << patch_grid.computeMinChLength() << endl;
cout << " dt  =  " << dt << endl;

//real dt = cfl_target*patch.dx()/(max(u_jet, u_far) + sqrt(PerfectGas::gamma()*PerfectGas::R()*max(T_jet,T)));

// Initialize
real init_var[5];
PerfectGas::primitiveToConservative(p, T, u_far, 0, 0, init_var);
patch_grid.setFieldToConst(0, init_var);

//cout << NI*NJ*NK << " cells" << endl;

JetFlux flux(u_jet, u_far, p_jet, p, T_jet, T, L, 2*L, L, 2*L, false);

//PerfectGas::primitiveToConservative(p, T, u_far, 1e-2*u_far, 2e-2*u_far, init_var);
//for (size_t i = 0; i < NI; ++i) {
//  for (size_t j = 0; j < NJ; ++j) {
//    for (size_t k = 0; k < NK; ++k) {
//      patch.setVar(0, i, j, k, init_var);
//    }
//  }
//}
