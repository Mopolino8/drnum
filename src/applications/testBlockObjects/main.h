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
// + enGrid is distributed in the hope that it will be useful,            +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with enGrid. If not, see <http://www.gnu.org/licenses/>.       +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef TESTBLOCKOBJECTS_H
#define TESTBLOCKOBJECTS_H

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

#ifdef GPU
#include "iterators/gpu_cartesianiterator.h"
#else
#include "iterators/cartesianiterator.h"
#endif

#include "rungekutta.h"
#include "iteratorfeeder.h"
#include "cubeincartisianpatch.h"

#include "blockobject.h"
#include "cartboxobject.h"
#include "sphereobject.h"
#include "cylinderobject.h"
#include "coneobject.h"
#include "combiobjectand.h"
#include "combiobjector.h"
#include "combiobjectandnot.h"


class EaFlux
{

protected:

  //typedef Upwind1 reconstruction_t;
  //typedef Upwind2<SecondOrder>                           reconstruction_t;

  typedef Upwind2<VanAlbada>                             reconstruction_t;
  typedef AusmPlus<reconstruction_t, PerfectGas>         euler_t;
  //typedef KNP<reconstruction_t, PerfectGas>              euler_t;
  typedef CompressibleSlipFlux<Upwind1, PerfectGas>      wall_t;
  typedef CompressibleViscFlux<PerfectGas>               viscous_t;
  typedef CompressibleFarfieldFlux<Upwind1, PerfectGas>  farfield_t;

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
  m_FarfieldFlux.yWallP(P, i, j, k, x, y, z, A, flux);
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
  m_FarfieldFlux.yWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename PATCH>
inline void EaFlux::zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarfieldFlux.zWallM(P, i, j, k, x, y, z, A, flux);
}


class EaWFluxZM : public EaFlux
{

public:

  EaWFluxZM(real u, real v, real p, real T) : EaFlux(u, v, p, T) {}
  EaWFluxZM() {}

  template <typename PATCH>
  CUDA_DH void zWallM(PATCH *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
  {
    m_WallFlux.zWallM(P, i, j, k, x, y, z, A, flux);
  }

};




//vector<CubeInCartisianPatch> setupCubes(PatchGrid patch_grid)
//{
//  vector<CubeInCartisianPatch> cubes;

//  // Building 1
//  {
//    vec3_t x1(-42.30365, -5, 0);
//    vec3_t x2(-13.46731, 9.07423, 14.29426);
//    {
//      //size_t i_patch = 16;
//      size_t i_patch = 60;
//      CartesianPatch* cart_patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_patch));
//      CubeInCartisianPatch cube(cart_patch);
//      cube.setRange(x1, x2);
//      cubes.push_back(cube);
//    }
//  }

//  // Bulding 2
//  {
//    vec3_t x1(-5, -5, 0);
//    vec3_t x2(5, 5, 65);
//    {
//      //size_t i_patch = 28;
//      size_t i_patch = 85;
//      CartesianPatch* cart_patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_patch));
//      CubeInCartisianPatch cube(cart_patch);
//      cube.setRange(x1, x2);
//      cubes.push_back(cube);
//    }
//    {
//      //size_t i_patch = 29;
//      size_t i_patch = 86;
//      CartesianPatch* cart_patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_patch));
//      CubeInCartisianPatch cube(cart_patch);
//      cube.setRange(x1, x2);
//      cubes.push_back(cube);
//    }
//  }

//  // Bulding 3
//  {
//    vec3_t x1(22.04841, -5, 0);
//    vec3_t x2(34.91528, 18.73290, 84.33866);
//    {
//      //size_t i_patch = 40;
//      size_t i_patch = 110;
//      CartesianPatch* cart_patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_patch));
//      CubeInCartisianPatch cube(cart_patch);
//      cube.setRange(x1, x2);
//      cubes.push_back(cube);
//    }
//    {
//      //size_t i_patch = 41;
//      size_t i_patch = 111;
//      CartesianPatch* cart_patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_patch));
//      CubeInCartisianPatch cube(cart_patch);
//      cube.setRange(x1, x2);
//      cubes.push_back(cube);
//    }
//    {
//      //size_t i_patch = 42;
//      size_t i_patch = 112;
//      CartesianPatch* cart_patch = dynamic_cast<CartesianPatch*>(patch_grid.getPatch(i_patch));
//      CubeInCartisianPatch cube(cart_patch);
//      cube.setRange(x1, x2);
//      cubes.push_back(cube);
//    }
//  }

//  return cubes;
//}


void run()
{

  real Ma             = 0.15;
  real p              = 1.0e5;
  real T              = 300;
  real uabs           = Ma*sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
  real alpha          = 10.0;
  real L              = 10.0;
  real time           = L/uabs;
  real cfl_target     = 1.0;
  real t_write        = 0;
  real write_interval = 1.0*time;
  real total_time     = 100*time;
  bool write          = true;
  bool write_blockobj = true;


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
  PerfectGas::primitiveToConservative(p, T, u, v, 0, init_var);
  patch_grid.setFieldToConst(0, init_var);

  // Define blockobject geometries here:

  // #define TEST_1
  // #define TEST_2
  // #define TEST_3
  // #define TEST_4
  #define TEST_5

#ifdef TEST_1

  // Test 1: single sphere
  SphereObject object;
  object.setParams(5., 5., 5., 3.);

#elif TEST_2

  // some spheres
  SphereObject obj_sph_1;
  SphereObject obj_sph_2;
  SphereObject obj_sph_3;
  obj_sph_1.setParams(4., 4., 4., 2.5);
  obj_sph_2.setParams(6., 6., 6., 2.5);
  obj_sph_3.setParams(4., 5., 6., 2.5);
  // some cylinders
  CylinderObject obj_cyl_1;
  obj_cyl_1.setParams(1., 5., 5.,
                      6., 0., 0.,
                      1.);

  // Combis of two
  // CombiObjectOr object(&obj_sph_1, &obj_sph_2);
  // CombiObjectAnd object(&obj_sph_1, &obj_sph_2);
  // CombiObjectAndNot object(&obj_sph_1, &obj_sph_2);

  // Combis of three
  // test 3.1
  // CombiObjectAnd object_inter1(&obj_sph_1, &obj_sph_2);
  // CombiObjectAndNot object(&obj_sph_3, &object_inter1);
  // test 3.2
  CombiObjectAnd object_inter1(&obj_sph_1, &obj_sph_2);
  CombiObjectAnd object_inter2(&obj_sph_2, &obj_sph_3);
  CombiObjectAnd object_inter3(&obj_sph_3, &obj_sph_1);
  CombiObjectOr object_inter4(&object_inter1, &object_inter2);
  CombiObjectOr object_inter5(&object_inter4, &object_inter3);
  CombiObjectAndNot object(&object_inter5, &obj_cyl_1);
#endif
#ifdef TEST_3

  // Model a desk
  CartboxObject plate;
  CartboxObject leg_1;
  CartboxObject leg_2;
  CartboxObject leg_3;
  CartboxObject leg_4;
  SphereObject vase_body;
  CylinderObject vase_inner;

  real x_p_min = 2.0;
  real x_p_max = 7.0;
  real y_p_min = 3.0;
  real y_p_max = 7.0;
  real leg_thick = 0.5;

  plate.setParams(x_p_min, x_p_max,
                  y_p_min, y_p_max,
                  2.5, 3.0);

  leg_1.setParams(x_p_min,             x_p_min + leg_thick,
                  y_p_min,             y_p_min + leg_thick,
                  -1.0, 2.75);

  leg_2.setParams(x_p_min,             x_p_min + leg_thick,
                  y_p_max - leg_thick, y_p_max,
                  -1.0, 2.75);

  leg_3.setParams(x_p_max - leg_thick, x_p_max,
                  y_p_min,             y_p_min + leg_thick,
                  -1.0, 2.75);

  leg_4.setParams(x_p_max - leg_thick, x_p_max,
                  y_p_max - leg_thick, y_p_max,
                  -1.0, 2.75);

  vase_body.setParams(4., 5., 3.8, 1.0);
  vase_inner.setParams(4., 5., 3.3,
                       0., 0., 100.,
                       0.3);

  // Combi objects
  // desk; plate and legs


  //  /// @todo How about a multiple boolean OR ?
  //  /// @todo Or a general class with boolean Ops
  //  CombiObjectOr h1(&leg_1, &leg_2);
  //  CombiObjectOr h2(&h1, &leg_3);
  //  CombiObjectOr h3(&h2, &leg_4);
  //  CombiObjectOr desk(&plate, &h3);
  ////  CombiObjectOr desk(&plate,
  ////        &CombiObjectOr(&leg_1,
  ////          &CombiObjectOr(&leg_2
  ////            &CombiObjectOr(&leg_3, &leg_4))));
  //  // vase
  //  CombiObjectAndNot vase(&vase_body, &vase_inner);
  //  // together
  //  CombiObjectOr object(&desk, &vase);


  CombiObjectOr desk (&plate);
  desk.includeObject(&leg_1);
  desk.includeObject(&leg_2);
  desk.includeObject(&leg_3);
  desk.includeObject(&leg_4);

  CombiObjectAndNot vase (&vase_body);
  vase.includeObject(&vase_inner);

  CombiObjectOr object(&desk, &vase);

#endif
#ifdef TEST_4

  ConeObject cone;
  cone.setParams(1., 5., 5.,
                 8., 0., 0.,
                 3., 1.);
  CylinderObject drill_1;
  drill_1.setParams(5., 5., -3.,
                    0., 0., 100.,
                    1.);
  CylinderObject drill_2;
  drill_2.setParams(4., -3., 5.,
                    0., 100., 0.,
                    1.);
  CombiObjectAndNot object (&cone);
  object.includeObject(&drill_1);
  object.includeObject(&drill_2);

#endif
#ifdef TEST_5

  CylinderObject object;
  object.setParams(2., 5., 5.,
                   6., 0., 0.,
                   1.);

#endif

      // Transform in block object on patch_grid
      BlockObject block_object(&patch_grid);
  block_object.setGreyResolution(5); // overwrite default
  block_object.update(&object);
  // Produce output of block object for visualization (abuse field 0 , variable 0)
  block_object.setLayerIndexToVar(0, 0);
  // Write plotteable output
  if (write_blockobj) {
    /// @todo Abuse field=0, var=0 to write file with blockobjects. Define other write method later.
    patch_grid.writeToVtk(0, "VTK/blockobjects", CompressibleVariables<PerfectGas>(), 0);
  }

  if (write) {
    patch_grid.writeToVtk(0, "VTK/step", CompressibleVariables<PerfectGas>(), 0);
  }

  exit(-1); /// meshing!!!

  EaFlux    flux_std(u, v, p, T);
  EaWFluxZM flux_wzm(u, v, p, T);

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

#ifdef GPU
  GPU_CartesianIterator<5, EaFlux>    iterator_std(flux_std);
  GPU_CartesianIterator<5, EaWFluxZM> iterator_wzm(flux_wzm);
#else
  CartesianIterator<5, EaFlux>    iterator_std(flux_std);
  CartesianIterator<5, EaWFluxZM> iterator_wzm(flux_wzm);
#endif
  iterator_std.setCodeString(CodeString("fx fy fz far far far far far  far 0"));
  iterator_wzm.setCodeString(CodeString("fx fy fz far far far far wall far 0"));

  IteratorFeeder iterator_feeder;
  iterator_feeder.addIterator(&iterator_std);
  iterator_feeder.addIterator(&iterator_wzm);
  iterator_feeder.feed(patch_grid);

  //  vector<CubeInCartisianPatch> cubes = setupCubes(patch_grid);
  //  for (size_t i = 0; i < cubes.size(); ++i) {
  //    runge_kutta.addPostOperation(&cubes[i]);
  //  }

  runge_kutta.addIterator(&iterator_std);
  runge_kutta.addIterator(&iterator_wzm);

  int write_counter = 0;
  int iter = 0;
  real t = 0;

  cout << "std:" << iterator_std.numPatches() << endl;
  cout << "wzm:" << iterator_wzm.numPatches() << endl;

#ifdef GPU
  iterator_std.updateDevice();
  iterator_wzm.updateDevice();
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
      iterator_wzm.updateHost();
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
      cout << iter << " iterations,  t=" << t/time << "*L/u_oo,  dt: " << dt;
      cout << "  CFL: " << CFL_max;
      cout << "  max: " << max_norm_allpatches << "  L2: " << l2_norm_allpatches;
      cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
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
  iterator_wzm.updateHost();
#endif

  if (write) {
    patch_grid.writeToVtk(0, "VTK/final", CompressibleVariables<PerfectGas>(), -1);
  }
}

#endif // TESTBLOCKOBJECTS_H
