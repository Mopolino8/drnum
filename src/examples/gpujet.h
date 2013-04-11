#ifndef GPUJET_H
#define GPUJET_H

#include "examples/jetflux.h"
#include "iterators/gpu_cartesianiterator.h"
#include "rungekutta.h"

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

void run()
{

  #include "jet_common.h"

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

  PatchGrid patch_grid;
  GPU_CartesianIterator<5, JetFlux> iterator(patch_grid, flux);
  iterator.addPatch(&patch);
  runge_kutta.addIterator(&iterator);

  int write_counter = 0;
  int iter = 0;
  real t = 0;
  write(patch, "testrun", write_counter);

  //cout << "Press <ENTER> to start!";
  //cin.get();

  startTiming();

  while (t < total_time) {

    runge_kutta(dt);
    t += dt;
    t_write += dt;

    real x = 0.5*patch.dx();
    real rho_min = 1000;
    real rho_max = 0;

    bool print_full = false;

    if (t_write >= write_interval) {
      print_full = true;
    }

    if (print_full) {
      real CFL_max = 0;
      iterator.updateHost();
      for (size_t i = 0; i < NI; ++i) {
        real y = 0.5*patch.dy();
        for (size_t j = 0; j < NJ; ++j) {
          real z = 0.5*patch.dz();
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
            z += patch.dz();
          }
          y += patch.dy();
        }
        x += patch.dx();
      }
      dt *= cfl_target/CFL_max;
    }
    if (t_write >= write_interval) {
      ++write_counter;
      iterator.updateHost();
      write(patch, "testrun", write_counter);
      t_write -= write_interval;
    }
    ++iter;
    cout << iter << " iterations,  t=" << t/time << "*L/u_jet,  dt: " << dt;;
    if (print_full) {
      real max_norm, l2_norm;
      patch.computeVariableDifference(0, 0, 1, 0, max_norm, l2_norm);
      cout << "  max: " << max_norm << "  L2: " << l2_norm;
      cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max;
    }
    cout << endl;
  }

  stopTiming();
  cout << iter << " iterations" << endl;

  iterator.updateHost();
  write(patch, "testrun_final", -1);
  //gpu_iterator.gpuInfo();
}

#endif // GPUJET_H
