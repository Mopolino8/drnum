#ifndef GPUJET_H
#define GPUJET_H

#include "examples/jetflux.h"
#include "iterators/gpu_cartesianiterator.h"
#include "rungekutta.h"
#include "iteratorfeeder.h"

void run()
{

  //string grid_file_name = "grid/jet_dual_2_mJ12_mK12.grid";
  //string grid_file_name = "grid/jet_dual_2.grid";
  string grid_file_name = "grid/bla.grid";

  #include "jet_pg_common.h"

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

  GPU_CartesianIterator<5, JetFlux> iterator(flux);
  iterator.setCodeString(CodeString("fx1 fy1 fz1 bx1 bX1 by1 bY1 bz1 bZ1 s1"));

  IteratorFeeder iterator_feeder;
  iterator_feeder.addIterator(&iterator);
  iterator_feeder.feed(patch_grid);

  runge_kutta.addIterator(&iterator);

  int write_counter = 0;
  int iter = 0;
  real t = 0;
  patch_grid.writeToVtk(0, "testrun", CompressibleVariables<PerfectGas>(), write_counter);
  iterator.updateDevice();

  //cout << "Press <ENTER> to start!";
  //cin.get();

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

      iterator.updateHost();

      for (size_t i_p = 0; i_p < patch_grid.getNumPatches(); i_p++) {
        // Patch* patch = patch_grid.getPatch(i_p);
        CartesianPatch& patch = *(dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_p)));
        size_t NI = patch.sizeI();
        size_t NJ = patch.sizeJ();
        size_t NK = patch.sizeK();

        for (size_t i = 0; i < NI; ++i) {
          for (size_t j = 0; j < NJ; ++j) {
            for (size_t k = 0; k < NK; ++k) {
              real p, u, v, w, T, var[5];
              patch.getVar(0, i, j, k, var);
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
      cout << iter << " iterations,  t=" << t/time << "*L/u_jet,  dt: " << dt;;
      cout << t/time << "  dt: " << dt << "  CFL: " << CFL_max;
      cout << "  max: " << max_norm_allpatches << "  L2: " << l2_norm_allpatches;
      cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
      ++write_counter;
      patch_grid.writeToVtk(0, "testrun", CompressibleVariables<PerfectGas>(), write_counter);
      t_write -= write_interval;
    } else {
      ++iter;
      cout << iter << " iterations,  t=" << t/time << "*L/u_jet,  dt: " << dt << endl;
    }
  }

  stopTiming();
  cout << iter << " iterations" << endl;

  iterator.updateHost();
  patch_grid.writeToVtk(0, "testrun_final", CompressibleVariables<PerfectGas>(), -1);
  //gpu_iterator.gpuInfo();
}

#endif // GPUJET_H
