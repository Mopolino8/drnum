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

#ifndef JETDEMO_MAIN_H
#define JETDEMO_MAIN_H

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/minmod.h"
#include "reconstruction/secondorder.h"
#include "fluxes/knp.h"
//#include "fluxes/vanleer.h"
#include "fluxes/ausmplus.h"
//#include "fluxes/ausmdv.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "fluxes/compressibleslipflux.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"
#include "externalexchangelist.h"

#ifdef GPU
#include "iterators/gpu_cartesianiterator.h"
#else
#include "iterators/cartesianiterator.h"
#endif

#include "rungekutta.h"
#include "iteratorfeeder.h"
#include "cubeincartisianpatch.h"

#include "configmap.h"

typedef Upwind2<5, VanAlbada> reconstruction_t;

#include "jetflux.h"
#include "eaflux.h"

void run()
{
  dim_t<5> dim;

  // control files
  ConfigMap config;
  config.addDirectory("control");

  real Ma_far         = config.getValue<real>("Ma-far");
  real Ma_jet         = config.getValue<real>("Ma-jet");
  real p_far          = config.getValue<real>("p-far");
  real p_jet          = config.getValue<real>("p-jet");
  real T_far          = config.getValue<real>("T-far");
  real T_jet          = config.getValue<real>("T-jet");
  real u_far          = Ma_far*sqrt(PerfectGas::gamma()*PerfectGas::R()*T_far);
  real u_jet          = Ma_jet*sqrt(PerfectGas::gamma()*PerfectGas::R()*T_jet);
  real L              = 2*config.getValue<real>("radius");
  real time           = L/u_jet;
  real cfl_target     = config.getValue<real>("CFL");
  real t_write        = 0;
  real write_interval = config.getValue<real>("write-interval")*time;
  real total_time     = config.getValue<real>("total-time")*time;
  int  jet_patch_id   = config.getValue<int>("jet-patch");
  bool mesh_preview   = config.getValue<bool>("mesh-preview");
  real scale          = config.getValue<real>("scale");
  int  thread_limit   = 0;

  if (config.exists("thread-limit")) {
    thread_limit = config.getValue<int>("thread-limit");
  }

#ifdef GPU
  int  cuda_device    = config.getValue<int>("cuda-device");
#endif

  // Patch grid
  PatchGrid patch_grid;
  //.. general settings (apply to all subsequent patches)
  patch_grid.setNumberOfFields(3);
  patch_grid.setNumberOfVariables(5);
  patch_grid.defineVectorVar(1);
  patch_grid.setInterpolateData();
  patch_grid.setNumSeekLayers(2);  /// @todo check default = 2
  patch_grid.setTransferType("padded_direct");
  patch_grid.readGrid("patches/standard.grid", scale);
  //patch_grid.readGrid("patches/V1");
  patch_grid.computeDependencies(true);

  // Time step
  real ch_speed1 = u_jet + sqrt(PerfectGas::gamma()*PerfectGas::R()*T_jet);
  real ch_speed2 = u_far + sqrt(PerfectGas::gamma()*PerfectGas::R()*T_far);
  real ch_speed  = max(ch_speed1, ch_speed2);
  real dt        = cfl_target*patch_grid.computeMinChLength()/ch_speed;

  cout << " patch_grid.computeMinChLength() = " << patch_grid.computeMinChLength() << endl;
  cout << " dt  =  " << dt << endl;

  // Initialize
  real init_var[5];
  PerfectGas::primitiveToConservative(p_far, T_far, 0, u_far, u_far, init_var);
  patch_grid.setFieldToConst(0, init_var);

  if (write) {
    patch_grid.writeToVtk(0, "VTK/step", CompressibleVariables<PerfectGas>(), 0);
  }

  if (mesh_preview) {
    exit(EXIT_SUCCESS);
  }

  JetFlux flux;
  flux.setup(u_jet, u_far, p_jet, p_far, T_jet, T_far, 0.5*L, jet_patch_id, max(Ma_far, Ma_jet) > 1.0);

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

#ifdef GPU
  GPU_CartesianIterator<5, JetFlux>  iterator_std(flux, cuda_device, thread_limit);
#else
  BUG;
  CartesianIterator<5, JetFlux>  iterator_std(flux);
#endif
  iterator_std.setCodeString(CodeString("fx fy fz far far far far far far 0"));

  IteratorFeeder iterator_feeder;
  iterator_feeder.addIterator(&iterator_std);
  iterator_feeder.feed(patch_grid);

  //vector<CubeInCartisianPatch> cubes = setupCubes(patch_grid);
  //for (size_t i = 0; i < cubes.size(); ++i) {
  //  runge_kutta.addPostOperation(&cubes[i]);
  //}
  //CartesianCyclicCopy<5> cyclic(&patch_grid);
  //runge_kutta.addPostOperation(&cyclic);

  runge_kutta.addIterator(&iterator_std);

  int write_counter = 0;
  int iter = 0;
  real t = 0;

#ifdef GPU
  iterator_std.updateDevice();
#endif

  startTiming();


  while (t < total_time) {

    runge_kutta(dt);
    t += dt;
    t_write += dt;

    if (t_write >= write_interval) {

      // Do some diagnose on patches
      real CFL_max = 0;
      real rho_min = 1000;
      real rho_max = 0;
      real max_norm_allpatches = 0.;
      real l2_norm_allpatches;
      real ql2_norm_allpatches = 0.;

#ifdef GPU
      runge_kutta.copyDonorData(0);
      iterator_std.updateHost();
#endif

      for (size_t i_p = 0; i_p < patch_grid.getNumPatches(); i_p++) {
        CartesianPatch& patch = *(dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_p)));
        size_t NI = patch.sizeI();
        size_t NJ = patch.sizeJ();
        size_t NK = patch.sizeK();

        for (size_t i = 0; i < NI; ++i) {
          for (size_t j = 0; j < NJ; ++j) {
            for (size_t k = 0; k < NK; ++k) {
              real p, u, v, w, T, var[5];
              patch.getVar(dim, 0, i, j, k, var);
              rho_min = min(var[0], rho_min);
              rho_max = max(var[0], rho_max);
              PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
              real a = sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
              CFL_max = max(CFL_max, fabs(u)*dt/patch.dx());
              CFL_max = max(CFL_max, fabs(u+a)*dt/patch.dx());
              CFL_max = max(CFL_max, fabs(u-a)*dt/patch.dx());
              countSqrts(1);
              countFlops(10);
            }
          }
        }
        real max_norm, l2_norm;
        patch.computeVariableDifference(0, 0, 1, 0, max_norm, l2_norm);
        if (max_norm > max_norm_allpatches) {
          max_norm_allpatches = max_norm;
        }
        ql2_norm_allpatches += l2_norm * l2_norm;
      }
      ql2_norm_allpatches /= patch_grid.getNumPatches();
      l2_norm_allpatches = sqrt(ql2_norm_allpatches);
      ++iter;
      cout << iter << " iterations,  t=" << t/time << "*L/u_oo,  dt: " << dt;
      cout << "  CFL: " << CFL_max;
      cout << "  max: " << max_norm_allpatches << "  L2: " << l2_norm_allpatches;
      cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
      printTiming();

      ConfigMap config;
      config.addDirectory("control");
      cfl_target     = config.getValue<real>("CFL");
      write_interval = config.getValue<real>("write-interval")*time;
      total_time     = config.getValue<real>("total-time")*time;

      dt *= cfl_target/CFL_max;

      ++write_counter;
      if (write) {
        patch_grid.writeToVtk(0, "VTK/step", CompressibleVariables<PerfectGas>(), write_counter);
      }
      t_write -= write_interval;
      if (t > 80.0) {
        write_interval = 0.05;
      }
    } else {
      ++iter;
      cout << iter << " iterations,  t=" << t << " = " << t/time << "*L/u_oo,  dt: " << dt << endl;
    }
  }

  stopTiming();
  cout << iter << " iterations" << endl;

#ifdef GPU
  runge_kutta.copyDonorData(0);
  iterator_std.updateHost();
#endif

  if (write) {
    patch_grid.writeToVtk(0, "VTK/final", CompressibleVariables<PerfectGas>(), -1);
  }
}

#endif // EXTERNAL_AERO_H
