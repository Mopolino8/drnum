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
#include "fluxes/kt.h"
#include "fluxes/knp.h"
#include "fluxes/vanleer.h"
#include "fluxes/ausmplus.h"
#include "fluxes/ausmdv.h"
#include "fluxes/roe.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressiblesmagorinskyflux.h"
#include "fluxes/compressibleslipflux.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"

#ifdef GPU
#include "cartesiancycliccopy.h"
#endif

#include "externalexchangelist.h"

#ifdef GPU
#include "iterators/gpu_cartesianiterator.h"
#else
#include "iterators/cartesianiterator.h"
#endif

#include "rungekutta.h"
#include "iteratorfeeder.h"
#include "cubeincartisianpatch.h"

#include <QTime>

#include "configmap.h"

template <typename TReconstruction>
class EaFlux
{

protected:

  //typedef AusmDV<5, reconstruction_t, PerfectGas> euler_t;
  //typedef VanLeer<5, TReconstruction, PerfectGas> euler_t;
  typedef Roe<5, TReconstruction, PerfectGas> euler_t;
  //typedef AusmPlus<5, TReconstruction, PerfectGas> euler_t;
  //typedef KT<5, 10000, TReconstruction, PerfectGas> euler_t;
  //typedef KNP<5, TReconstruction, PerfectGas> euler_t;

  typedef CompressibleSlipFlux<5, Upwind1<5>, PerfectGas>     wall_t;
  typedef CompressibleSmagorinskyFlux<5, 2000, PerfectGas>    viscous_t;
  typedef CompressibleFarfieldFlux<5, Upwind1<5>, PerfectGas> farfield_t;
  typedef CompressibleFlux<5, PerfectGas>                     split_t;

  TReconstruction m_Reconstruction;
  euler_t          m_EulerFlux;
  viscous_t        m_ViscFlux;
  farfield_t       m_FarfieldFlux;
  wall_t           m_WallFlux;
  split_t          m_SplitFlux;
  bool             m_Inviscid;


public: // methods

  EaFlux(real u, real v, real p, real T, bool inviscid);
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

  template <typename PATCH> CUDA_DH void splitFlux(PATCH *P, splitface_t sf, real* flux);

};


template <typename TReconstruction>
EaFlux<TReconstruction>::EaFlux(real u, real v, real p, real T, bool inviscid)
{
  m_FarfieldFlux.setFarfield(p, T, u, v, 0);
  m_Inviscid = inviscid;
}

template <typename TReconstruction>
EaFlux<TReconstruction>::EaFlux()
{
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::xField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.xField(patch, i, j, k, x, y, z, A, flux);
  if (!m_Inviscid) m_ViscFlux.xField(patch, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::yField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.yField(patch, i, j, k, x, y, z, A, flux);
  if (!m_Inviscid) m_ViscFlux.yField(patch, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::zField
(
  PATCH *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux.zField(patch, i, j, k, x, y, z, A, flux);
  if (!m_Inviscid) m_ViscFlux.zField(patch, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::xWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.xWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::yWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarfieldFlux.yWallP(P, i, j, k, x, y, z, A, flux);
  m_WallFlux.yWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::zWallP(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::xWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.xWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::yWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarfieldFlux.yWallM(P, i, j, k, x, y, z, A, flux);
  m_WallFlux.yWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename TReconstruction>
template <typename PATCH>
inline void EaFlux<TReconstruction>::splitFlux(PATCH *P, splitface_t sf, real* flux)
{
  m_SplitFlux.splitFlux(P, sf, flux);
  if (!m_Inviscid) m_ViscFlux.splitFlux(P, sf, flux);
}

void sync(Patch *patch, ExternalExchangeList *of2dn_list, ExternalExchangeList *dn2of_list, Barrier *barrier, SharedMemory *shmem, bool &write_flag, bool &stop_flag, real &dt)
{
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
  barrier->wait();
  dn2of_list->ipcSend();
  barrier->wait();
  int write = 0;
  int stop = 0;
  shmem->readValue("write", &write);
  shmem->readValue("stop", &stop);
  shmem->readValue("dt", &dt);
  if (write) {
    write_flag = true;
  } else {
    write_flag = false;
  }
  if (stop) {
    stop_flag = true;
    cout << "Received stop signal from client." << endl;
  } else {
    stop_flag = false;
  }
  of2dn_list->ipcReceive();

  for (int i = 0; i < of2dn_list->size(); ++i) {
    real var[5];
    p = of2dn_list->data(0, i);
    u = of2dn_list->data(1, i);
    v = of2dn_list->data(2, i);
    w = of2dn_list->data(3, i);
    T = of2dn_list->data(4, i);
    PerfectGas::primitiveToConservative(p, T, u, v, w, var);
    patch->setVar(dim, 0, of2dn_list->index(i), var);
  }

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
  real time           = L/uabs;
  real cfl_target     = config.getValue<real>("CFL-number");
  real t_write        = 0;
  real write_interval = config.getValue<real>("write-interval")*time;
  real total_time     = config.getValue<real>("total-time")*time;
  bool mesh_preview   = config.getValue<bool>("mesh-preview");
  bool code_coupling  = config.getValue<bool>("code-coupling");
  real scale          = config.getValue<real>("scale");

  int thread_limit   = 0;
  if (config.exists("thread-limit")) {
    thread_limit = config.getValue<int>("thread-limit");
  }
#ifdef GPU
  int  cuda_device    = config.getValue<int>("cuda-device");
#endif

  bool write_flag     = false;
  bool stop_flag      = false;

  bool start_from_zero = false;
  if (config.exists("start-from-zero")) {
    start_from_zero = config.getValue<bool>("start-from-zero");
  }

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
  patch_grid.readGrid("patches/standard.grid", scale);
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
    patch_grid.writeToVtk(0, "VTK-drnum/step", CompressibleVariables<PerfectGas>(), 0);
  }

  if (mesh_preview) {
    exit(EXIT_SUCCESS);
  }


  RungeKutta runge_kutta;
  {
    int num_rk_steps = config.getValue<int>("num-rk-steps");
    list<real> alpha;
    real coeff = 1.0;
    for (int i = 0; i < num_rk_steps; ++i) {
      alpha.push_front(coeff);
      coeff *= 0.5;
    }
    for (list<real>::iterator i = alpha.begin(); i != alpha.end(); ++i) {
      runge_kutta.addAlpha(*i);
    }
  }

  QString reconstruction = config.getValue<QString>("reconstruction");

  PatchIterator *iterator;

  bool inviscid = config.getValue<bool>("inviscid");

  if (reconstruction == "second-order") {
    EaFlux<Upwind2<5, SecondOrder> > flux(u, v, p, T, inviscid);
#ifdef GPU
    iterator = new GPU_CartesianIterator<5, EaFlux<Upwind2<5, SecondOrder> > >(flux, cuda_device, thread_limit);
#else
    iterator = new CartesianIterator<5, EaFlux<Upwind2<SecondOrder> > >(flux);
#endif

  } else if (reconstruction == "minmod") {
    EaFlux<Upwind2<5, MinMod> > flux(u, v, p, T, inviscid);
#ifdef GPU
    iterator = new GPU_CartesianIterator<5, EaFlux<Upwind2<5, MinMod> > >(flux, cuda_device, thread_limit);
#else
    iterator = new CartesianIterator<5, EaFlux<Upwind2<MinMod> > >(flux);
#endif

  } else {
    EaFlux<Upwind1<5> > flux(u, v, p, T, inviscid);
#ifdef GPU
    iterator = new GPU_CartesianIterator<5, EaFlux<Upwind1<5> > >(flux, cuda_device, thread_limit);
#else
    iterator = new CartesianIterator<5, EaFlux<Upwind1<5> > >(flux);
#endif
  }

  iterator->setCodeString(CodeString("fx fy fz far far far far far  far 0"));

  IteratorFeeder iterator_feeder;
  iterator_feeder.addIterator(iterator);
  iterator_feeder.feed(patch_grid);

  runge_kutta.addIterator(iterator);

  int write_counter = 0;
  int iter = 0;
  real t = 0;


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
    sync(coupling_patch, of2dn_list, dn2of_list, barrier, shmem, write_flag, stop_flag, dt);
    write_interval = MAX_REAL;
    total_time = MAX_REAL;
    iterator->deactivatePatch(coupling_patch_id);
  }

  QString restart_file = config.getValue<QString>("restart-file");
  if (restart_file.toLower() != "none") {
    write_counter = restart_file.right(6).toInt();
    t = patch_grid.readData(0, "data/" + restart_file);
    if (start_from_zero) {
      t = 0;
    }
  }

  // set inside of bodies at rest (if requested)
  bool inside_at_rest = config.getValue<bool>("inside-at-rest");
  if (inside_at_rest) {
    for (size_t i_patch = 0; i_patch < patch_grid.getNumPatches(); ++i_patch) {
      Patch *P = patch_grid.getPatch(i_patch);
      for (size_t idx = 0; idx < P->variableSize(); ++idx) {
        if (P->isInsideCell(idx)) {
          P->getVariable(0, 1)[idx] = 0;
          P->getVariable(0, 2)[idx] = 0;
          P->getVariable(0, 3)[idx] = 0;
        }
      }
    }
  }

#ifdef GPU
  iterator->updateDevice();
#endif

  startTiming();

  while (t < total_time && !stop_flag) {

    QTime step_start = QTime::currentTime();
    runge_kutta(dt);
    int msecs_drnum = step_start.msecsTo(QTime::currentTime());
    real dt_new = dt;
    if (coupling_patch) {
      sync(coupling_patch, of2dn_list, dn2of_list, barrier, shmem, write_flag, stop_flag, dt_new);
    }
    int msecs_total = step_start.msecsTo(QTime::currentTime());
    real drnum_fraction = real(msecs_drnum)/real(msecs_total);
    t += dt;
    t_write += dt;

    if (t_write >= write_interval || write_flag) {

      // Do some diagnose on patches
      real CFL_max = 0;
      real rho_min = 1000;
      real rho_max = 0;
      real max_norm_allpatches = 0.;
      real l2_norm_allpatches;
      real ql2_norm_allpatches = 0.;

#ifdef GPU
      runge_kutta.copyDonorData(0);
      iterator->updateHost();
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
      cfl_target     = config.getValue<real>("CFL-number");
      write_interval = config.getValue<real>("write-interval")*time;
      total_time     = config.getValue<real>("total-time")*time;
      dt *= cfl_target/CFL_max;

      ++write_counter;
      if (write) {
        if (config.getValue<bool>("file-output")) {
          patch_grid.writeToVtk(0, "VTK-drnum/step", CompressibleVariables<PerfectGas>(), write_counter);
          patch_grid.writeData(0, "data/step", t, write_counter);
        }
      }
      t_write -= write_interval;
    } else {
      ++iter;
      cout << iter << " iterations,  t=" << t << ",  t=" << t/time << "*L/u_oo,  dt: " << dt << ",  DrNUM % of run-time: " << 100*drnum_fraction << endl;
    }
    if (coupling_patch) {
      dt = dt_new;
    }
  }

  stopTiming();
  cout << iter << " iterations" << endl;

#ifdef GPU
  runge_kutta.copyDonorData(0);
  iterator->updateHost();
#endif

  if (write) {
    patch_grid.writeToVtk(0, "VTK-drnum/final", CompressibleVariables<PerfectGas>(), -1);
  }
}

#endif // EXTERNAL_AERO_H
