#ifndef CPUJET_PG_H
#define CPUJET_PG_H

#include "examples/jetflux.h"
#include "iterators/cartesianiterator.h"
#include "iteratorfeeder.h"
#include "rungekuttapg1.h"


/// @todo write functions should go to patch_grid and patch or derived classes
void write(CartesianPatch &patch, QString file_name, int count)
{
  if (count >= 0) {
    QString num;
    num.setNum(count);
    while (num.size() < 6) {
      num = "0" + num;
    }
    file_name += "_" + num;
  }
  cout << "writing to file \"" << qPrintable(file_name) << ".vtr\"" << endl;
  CompressibleVariables<PerfectGas> cvars;
  patch.writeToVtk(0, file_name, cvars);
}

/// @todo write functions should go to patch_grid and patch or derived classes
void writePG(PatchGrid &patch_grid, QString file_name, int count)
{
  QString str_patch_index;
  QString str_patch_filename;
  // loop for patches iof patch_grid
  for (size_t i_p = 0; i_p < patch_grid.getNumPatches(); i_p++) {
    str_patch_index.setNum(i_p);
    while (str_patch_index.size() < 6) {
      str_patch_index = "0" + str_patch_index;
    }
    str_patch_filename = file_name + "_ip" + str_patch_index;
    CartesianPatch* patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_p));
    if(patch != NULL) {
      write(*patch, str_patch_filename, count);
    } else {
      cout << "Cant write data for other than CartesianPatch patches at present" << endl;
    }
  }
}


void run()
{
  //string grid_file_name = "grid/jet_single.grid";
  string grid_file_name = "grid/jet_dual_1.grid";

#include "jet_pg_common.h"

  RungeKuttaPG1 runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

  //CartesianStandardPatchOperation<5, JetFlux> operation(&patch, &flux);
  //CartesianStandardIterator iterator(&operation);
  //runge_kutta.addIterator(&iterator);

  CartesianIterator<5, JetFlux> iterator(flux);
  iterator.setCodeString(CodeString("fx1 fy1 fz1 bx1 bX1 by1 bY1 bz1 bZ1 s1"));

  IteratorFeeder iterator_feeder;
  iterator_feeder.addIterator(&iterator);
  iterator_feeder.feed(patch_grid);
  // iterator.addPatch(&patch);

  runge_kutta.addIterator(&iterator);

  int count = 0;
  int iter = 0;
  real t = 0;

  //CartesianPatch& patch_0 = *(dynamic_cast<CartesianPatch*>(patch_grid.getPatch(0)));
  writePG(patch_grid, "testrun", count);

  cout << "Press <ENTER> to start!";
  cin.get();

  startTiming();

  // Time loop
  while (t < total_time) {

    // Compute time step
    runge_kutta(dt);
    t += dt;
    t_write += dt;

    // Do some diagnose on patches
    real CFL_max = 0;
    real rho_min = 1000;
    real rho_max = 0;
    real max_norm_allpatches = 0.;
    real l2_norm_allpatches;
    real ql2_norm_allpatches = 0.;

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
    cout << t/time << "  dt: " << dt << "  CFL: " << CFL_max;
    cout << "  max: " << max_norm_allpatches << "  L2: " << l2_norm_allpatches;
    cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;

    // Write to output, if given time
    if (t_write >= write_interval) {
      ++count;
      dt *= cfl_target/CFL_max;
      writePG(patch_grid, "testrun", count);
      t_write -= write_interval;
    }

    // Count iterations
    ++iter;
  }

  stopTiming();
  cout << iter << " iterations" << endl;
}

#endif // CPUJET_PG_H
