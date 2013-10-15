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
#ifndef EXTERNAL_AERO_H
#define EXTERNAL_AERO_H

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
#include "cartesiancycliccopy.h"
#include "externalexchangelist.h"

#ifdef GPU
#include "iterators/gpu_cartesianiterator.h"
#else
#include "iterators/cartesianiterator.h"
#endif

#include "rungekutta.h"
#include "iteratorfeeder.h"
#include "cubeincartisianpatch.h"

// BEGIN testing
#include <QVector>
// END testing

#include "configmap.h"

class EaFlux
{

protected:

  //typedef Upwind1 reconstruction_t;
  typedef Upwind2<5, SecondOrder>                           reconstruction_t;

  //typedef Upwind2<5, VanAlbada>                             reconstruction_t;
  typedef AusmPlus<reconstruction_t, PerfectGas>         euler_t;
  //typedef KNP<reconstruction_t, PerfectGas>              euler_t;
  typedef CompressibleSlipFlux<5, Upwind1, PerfectGas>     wall_t;
  typedef CompressibleViscFlux<5, PerfectGas>               viscous_t;
  typedef CompressibleFarfieldFlux<5, Upwind1, PerfectGas>  farfield_t;

  reconstruction_t m_Reconstruction;
  euler_t          m_EulerFlux;
  viscous_t        m_ViscFlux;
  farfield_t       m_FarfieldFlux;
  wall_t           m_WallFlux;


public: // methods

  EaFlux(real u, real v, real p, real T);
  EaFlux();

  template <typename PATCH> CUDA_DH void xField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zField(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

  template <typename PATCH> CUDA_DH void xWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  template <typename PATCH> CUDA_DH void zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

};


EaFlux::EaFlux(real u, real v, real p, real T)
{
  m_FarfieldFlux.setFarfield(p, T, u, v, 0);
}

EaFlux::EaFlux()
{
}

template <typename PATCH>
inline void EaFlux::xField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.xField(patch, i, j, k, x, y, z, A, flux);
  m_ViscFlux.xField(patch, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::yField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.yField(patch, i, j, k, x, y, z, A, flux);
  m_ViscFlux.yField(patch, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::zField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.zField(patch, i, j, k, x, y, z, A, flux);
  m_ViscFlux.zField(patch, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::xWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.xWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarfieldFlux.yWallP(P, i, j, k, x, y, z, A, flux);
  m_WallFlux.yWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.xWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarfieldFlux.yWallM(P, i, j, k, x, y, z, A, flux);
  m_WallFlux.yWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallM(P, i, j, k, x, y, z, A, flux);
}



vector<CubeInCartisianPatch> setupCubes(PatchGrid patch_grid)
{
  vector<CubeInCartisianPatch> cubes;

  // plate
  {
    real L = 0.000305298288755695;
    real B = 0.1*L;
    real e = 0.005*L;
    vec3_t x1(-0.5*B, -0.5*L + e, -10*L);
    vec3_t x2( 0.5*B,  0.5*L + e,  10*L);
    {
      size_t i_patch = 12;
      CartesianPatch* cart_patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_patch));
      CubeInCartisianPatch cube(cart_patch);
      cube.setRange(x1, x2);
      cubes.push_back(cube);
    }
  }
  return cubes;
}

void sync(Patch *patch, ExternalExchangeList *of2dn_list, ExternalExchangeList *dn2of_list, Barrier *barrier)
{
  cout << "A" << endl;
  barrier->wait();

  dim_t<5> dim;
  patch->copyFieldToHost(0);
  real var[dim()];
  real p, T, u, v, w;
  for (int i = 0; i < dn2of_list->size(); ++i) {
    patch->getVar(dim, 0, dn2of_list->index(i), var);
    PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
    dn2of_list->data(0, i) = p;
    dn2of_list->data(1, i) = u;
    dn2of_list->data(2, i) = v;
    dn2of_list->data(3, i) = w;
    dn2of_list->data(4, i) = T;
  }
  dn2of_list->ipcSend();
  cout << "B" << endl;
  barrier->wait();
  cout << "C" << endl;
  of2dn_list->ipcReceive();

  //PerfectGas::primitiveToConservative(p, T, 0, 0, 0, var1);


  patch->copyFieldToDevice(0);
}

void run()
{
  dim_t<5> dim;

  // control files
  ConfigMap config;
  config.addDirectory("control");

  real Ma             = config.getValue<real>("Mach-number");
  real p              = config.getValue<real>("pressure");
  real T              = config.getValue<real>("temperature");
  real uabs           = Ma*sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
  real alpha          = 0.0;
  real L              = config.getValue<real>("reference-length");
  real time           = 1e-4/uabs;
  real cfl_target     = config.getValue<real>("CFL-number");
  real t_write        = 0;
  real write_interval = config.getValue<real>("write-interval")*time;
  real total_time     = config.getValue<real>("total-time")*time;
  bool mesh_preview   = config.getValue<bool>("mesh-preview");
  bool code_coupling  = config.getValue<bool>("code-coupling");

  alpha = M_PI*alpha/180.0;
  real u = uabs*cos(alpha);
  real v = uabs*sin(alpha);

  // Patch grid
  PatchGrid patch_grid;
  //.. general settings (apply to all subsequent patches)
  patch_grid.setNumberOfFields(3);
  patch_grid.setNumberOfVariables(5);
  patch_grid.defineVectorVar(1);
  patch_grid.setInterpolateData();
  patch_grid.setNumSeekLayers(2);  /// @todo check default = 2
  patch_grid.setTransferType("padded_direct");
  patch_grid.readGrid("patches/standard.grid");
  //patch_grid.readGrid("patches/V1");
  patch_grid.computeDependencies(true);

  // Time step
  real ch_speed = u + sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
  real dt       = cfl_target*patch_grid.computeMinChLength()/ch_speed;

  cout << " patch_grid.computeMinChLength() = " << patch_grid.computeMinChLength() << endl;
  cout << " dt  =  " << dt << endl;

  // Initialize
  real init_var[5];
  PerfectGas::primitiveToConservative(p, T, u, v, 0.01*u, init_var);
  patch_grid.setFieldToConst(0, init_var);

  if (write) {
    patch_grid.writeToVtk(0, "VTK/step", CompressibleVariables<PerfectGas>(), 0);
  }

  if (mesh_preview) {
    exit(EXIT_SUCCESS);
  }

  EaFlux    flux_std(u, v, p, T);

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

#ifdef GPU
  GPU_CartesianIterator<5, EaFlux>    iterator_std(flux_std);
#else
  CartesianIterator<5, EaFlux>    iterator_std(flux_std);
#endif
  iterator_std.setCodeString(CodeString("fx fy fz far far far far far  far 0"));

  IteratorFeeder iterator_feeder;
  iterator_feeder.addIterator(&iterator_std);
  iterator_feeder.feed(patch_grid);

  vector<CubeInCartisianPatch> cubes = setupCubes(patch_grid);
  for (size_t i = 0; i < cubes.size(); ++i) {
    runge_kutta.addPostOperation(&cubes[i]);
  }
  CartesianCyclicCopy<5> cyclic(&patch_grid);
  runge_kutta.addPostOperation(&cyclic);

  runge_kutta.addIterator(&iterator_std);

  int write_counter = 0;
  int iter = 0;
  real t = 0;

  cout << "std:" << iterator_std.numPatches() << endl;

  SharedMemory         *shmem = NULL;
  Barrier              *barrier = NULL;
  ExternalExchangeList *of2dn_list = NULL;
  ExternalExchangeList *dn2of_list = NULL;
  int                   coupling_patch_id = -1;
  Patch                *coupling_patch = NULL;

  if (code_coupling) {
    try {
      shmem = new SharedMemory(1, 32*1024*1024, true);
      barrier = new Barrier(2, true);
    } catch (IpcException E) {
      E.print();
    }
    of2dn_list = new ExternalExchangeList("of2dn", 5, NULL, shmem, barrier);
    dn2of_list = new ExternalExchangeList("dn2of", 5, NULL, shmem, barrier);
    int client_ready = 0;
    shmem->writeValue("client-ready", &client_ready);
    cout << "External code coupling has been enabled." << endl;
    cout << "waiting for client to connect ..." << endl;
    while (!client_ready) {
      shmem->readValue("client-ready", &client_ready);
    }
    barrier->wait();
    cout << "The client connection has been established." << endl;
    barrier->wait();
    of2dn_list->ipcReceive();
    dn2of_list->ipcReceive();
    barrier->wait();
    coupling_patch_id = config.getValue<int>("coupling-patch");
    of2dn_list->finalise(&patch_grid, coupling_patch_id);
    dn2of_list->finalise(&patch_grid, coupling_patch_id);
    coupling_patch = patch_grid.getPatch(coupling_patch_id);
    sync(coupling_patch, of2dn_list, dn2of_list, barrier);
    exit(-1);
  }

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
      cout << iter << " iterations,  t=" << t/time << "*L/u_oo,  dt: " << dt << endl;
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
