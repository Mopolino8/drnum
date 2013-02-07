#ifndef TWO_PATCHES_1_H
#define TWO_PATCHES_1_H

#include "shapes/halfspace.h"
#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/immersedboundaryreconstruction.h"
#include "fluxes/knp.h"
#include "fluxes/compressiblewallflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "boundary_conditions/compressibleeulerwall.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "iterators/cartesianstandarditerator.h"

#include "patchgrid.h"
#include "rungekuttapg1.h"

class MyFlux
{

  typedef ImmersedBoundaryReconstruction<Upwind2<VanAlbada>, HalfSpace, CompressibleEulerWall> reconstruction2_t;
  typedef ImmersedBoundaryReconstruction<Upwind1, HalfSpace, CompressibleEulerWall>            reconstruction1_t;
  typedef KNP<reconstruction1_t, PerfectGas> euler1_t;
  typedef KNP<reconstruction2_t, PerfectGas> euler2_t;
  typedef CompressibleWallFlux<reconstruction1_t, PerfectGas> wall_t;
  typedef CompressibleFarfieldFlux<reconstruction1_t, PerfectGas> farfield_t;

  reconstruction1_t*    m_Reconstruction1;
  reconstruction2_t*    m_Reconstruction2;
  euler1_t*             m_EulerFlux1;
  euler2_t*             m_EulerFlux2;
  wall_t*               m_WallFlux;
  farfield_t*           m_FarFlux;
  HalfSpace*            m_Shape;
  bool                  m_SecondOrder;


public: // methods

  MyFlux(real p, real T, real u, HalfSpace* shape);
  bool isInside(size_t i, size_t j, size_t k);

  void xField(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void yField(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void zField(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

  void xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);
  void zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux);

  void secondOrderOn()  { m_SecondOrder = true; }
  void secondOrderOff() { m_SecondOrder = false; }

};


MyFlux::MyFlux(real p, real T, real u, HalfSpace *shape)
{
  m_Shape = shape;
  m_Reconstruction1 = new reconstruction1_t(m_Shape);
  m_Reconstruction2 = new reconstruction2_t(m_Shape);
  m_EulerFlux1 = new euler1_t(m_Reconstruction1);
  m_EulerFlux2 = new euler2_t(m_Reconstruction2);
  m_FarFlux = new farfield_t(m_Reconstruction1);
  m_FarFlux->setFarfield(p, T, u, 0, 0);
  m_WallFlux = new wall_t(m_Reconstruction1);
  m_SecondOrder = false;
}

inline void MyFlux::xField
(
  CartesianPatch *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
  )
{
  if (m_SecondOrder) m_EulerFlux2->xField(patch, i, j, k, x, y, z, A, flux);
  else               m_EulerFlux1->xField(patch, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::yField
(
  CartesianPatch *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
  )
{
  if (m_SecondOrder) m_EulerFlux2->yField(patch, i, j, k, x, y, z, A, flux);
  else               m_EulerFlux1->yField(patch, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::zField
(
  CartesianPatch *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
  )
{
  if (m_SecondOrder) m_EulerFlux2->zField(patch, i, j, k, x, y, z, A, flux);
  else               m_EulerFlux1->zField(patch, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarFlux->xWallP(P, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux->yWallP(P, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux->zWallP(P, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarFlux->xWallM(P, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux->yWallM(P, i, j, k, x, y, z, A, flux);
}

inline void MyFlux::zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_WallFlux->zWallM(P, i, j, k, x, y, z, A, flux);
}

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
  // Test settings equivalent to wedge test wedge.h
  // Patch grid
  PatchGrid patch_grid;
  patch_grid.setInterpolateData();
  patch_grid.setTransferPadded();
  CartesianPatch* patches;
  patches = new CartesianPatch[2];  // this is just a test for 2 patches
  //.. patches 0 and 1
  real x_min, y_min, z_min, x_max, y_max, z_max;
  y_min = -0.01;
  y_max = 0.01;
  z_min = 0.;
  z_max = 1.;
  size_t N = 100;
  size_t NI = 3*N;
  size_t NJ = 1;
  size_t NK = 1*N;
  //.. patch 0
  // patches[0].setupAligned(0, -0.01, 0, 3, 0.01, 1);
  x_min = 0.;
  x_max = 1.5 + 0.02;  // check this: it should make two exact overlap layers
  size_t NI_0 = NI/2 + 2;
  patches[0].setNumberOfFields(2);
  patches[0].setNumberOfVariables(5);
  patches[0].setupAligned(x_min, y_min, z_min, x_max, y_max, z_max);
  patches[0].resize(NI_0, NJ, NK);
  patch_grid.insertPatch(&patches[0]);
  //.. patch 1
  // patches[0].setupAligned(0, -0.01, 0, 3, 0.01, 1);
  x_min = 1.5 - 0.02;
  x_max = 3.;
  size_t NI_1 = NI/2 + 2;
  patches[1].setNumberOfFields(2);
  patches[1].setNumberOfVariables(5);
  patches[1].setupAligned(x_min, y_min, z_min, x_max, y_max, z_max);
  patches[1].resize(NI_1, NJ, NK);
  patch_grid.insertPatch(&patches[1]);
  // Compute dependencies
  patch_grid.computeDependencies(true);
  //
  // Initialize
  real init_var[5];
  real zero_var[5];
  real Ma = 3;
  real p  = 1e5;
  real T  = 300;
  real u  = Ma*sqrt(PerfectGas::gamma()*PerfectGas::R()*T);

  PerfectGas::primitiveToConservative(p, T, u, 0, 0, init_var);
  PerfectGas::primitiveToConservative(p, T, 0, 0, 0, zero_var);
  for (size_t i = 0; i < NI_0; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        patches[0].setVar(0, i, j, k, init_var);
      }
    }
  }
  for (size_t i = 0; i < NI_1; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        patches[1].setVar(0, i, j, k, init_var);
      }
    }
  }

  //
  // Geometry
  ///@todo individual shape defs per patch not smart. Currently needed due to transformations: Reverse dependencies!
  typedef HalfSpace shape_t;
  shape_t* shapes;
  shapes = new shape_t[2];  // this is just a test for 2 patches
  shapes[0].setPoint(1,0,0);
  shapes[0].setNormal(-0.258819, 0, 0.9659258);
  shapes[0].transform(patches[0].getTransformation());
  shapes[1].setPoint(1,0,0);
  shapes[1].setNormal(-0.258819, 0, 0.9659258);
  shapes[1].transform(patches[1].getTransformation());
  //
  // Build up solver environment
  RungeKuttaPG1 runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

  MyFlux flux_0(p, T, u, &shapes[0]);  /// @todo ugly, even for testing
  MyFlux flux_1(p, T, u, &shapes[1]);  /// @todo ugly, even for testing
  CartesianStandardPatchOperation<5, MyFlux> operation_0(&patches[0], &flux_0);
  CartesianStandardPatchOperation<5, MyFlux> operation_1(&patches[1], &flux_1);
  CartesianStandardIterator iterator_0(&operation_0);
  CartesianStandardIterator iterator_1(&operation_1);
  iterator_0.setI1(1);
  iterator_1.setI1(1);

  runge_kutta.addIterator(&iterator_0);
  runge_kutta.addIterator(&iterator_1);

  runge_kutta.setPatchGrid(&patch_grid);
  runge_kutta.addSyncField(1);    /** @todo New field required, "1" OK? */

  real time           = 3.0/u;
  real dt             = 0.6*time/(3*N);
  real dt2            = dt;
  real t_switch       = 0*time;
  real t_write        = 0;
  real write_interval = time/2;
  real total_time     = 4*time;

  int count = 0;
  int iter = 0;
  real t = 0;

  write(patches[0], "testrun_0", count);
  write(patches[1], "testrun_1", count);

  cout << "Press <ENTER> to start!";
  cin.get();

  startTiming();

  while (t < total_time) {
    if (t > t_switch) {
      flux_0.secondOrderOn();
      flux_1.secondOrderOn();
      dt = dt2;
    }
    runge_kutta(dt);


    // Max CFL and max / min densities
    real CFL_max = 0;
    real rho_min = 1000;
    real rho_max = 0;

    // patch 0
    {
      real x = 0.5*patches[0].dx();
      for (size_t i = 0; i < NI; ++i) {
        real y = 0.5*patches[0].dy();
        for (size_t j = 0; j < NJ; ++j) {
          real z = 0.5*patches[0].dz();
          for (size_t k = 0; k < NK; ++k) {
            if (!shapes[0].isInside(x, y, z)) {
              real p, u, v, w, T, var[5];
              patches[0].getVar(0, i, j, k, var);
              rho_min = min(var[0], rho_min);
              rho_max = max(var[0], rho_max);
              PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
              real a = sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
              CFL_max = max(CFL_max, fabs(u)*dt/patches[0].dx());
              CFL_max = max(CFL_max, fabs(u+a)*dt/patches[0].dx());
              CFL_max = max(CFL_max, fabs(u-a)*dt/patches[0].dx());
              countSqrts(1);
              countFlops(10);
            } else {
              patches[0].setVar(0, i, j, k, zero_var);
            }
            z += patches[0].dz();
          }
          y += patches[0].dy();
        }
        x += patches[0].dx();
      }
    }

    // patch 1
    {
      real x = 0.5*patches[1].dx();
      for (size_t i = 0; i < NI; ++i) {
        real y = 0.5*patches[1].dy();
        for (size_t j = 0; j < NJ; ++j) {
          real z = 0.5*patches[1].dz();
          for (size_t k = 0; k < NK; ++k) {
            if (!shapes[1].isInside(x, y, z)) {
              real p, u, v, w, T, var[5];
              patches[1].getVar(0, i, j, k, var);
              rho_min = min(var[0], rho_min);
              rho_max = max(var[0], rho_max);
              PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
              real a = sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
              CFL_max = max(CFL_max, fabs(u)*dt/patches[1].dx());
              CFL_max = max(CFL_max, fabs(u+a)*dt/patches[1].dx());
              CFL_max = max(CFL_max, fabs(u-a)*dt/patches[1].dx());
              countSqrts(1);
              countFlops(10);
            } else {
              patches[1].setVar(0, i, j, k, zero_var);
            }
            z += patches[1].dz();
          }
          y += patches[1].dy();
        }
        x += patches[1].dx();
      }
    }

    // until here

    t += dt;
    t_write += dt;
    if (t_write >= write_interval) {
      ++count;
      write(patches[0], "testrun_0", count);
      write(patches[1], "testrun_1", count);
      t_write -= write_interval;
    }
    real max_norm_0, l2_norm_0;
    real max_norm_1, l2_norm_1;
    real max_norm, l2_norm;
    patches[0].computeVariableDifference(0, 0, 1, 0, max_norm_0, l2_norm_0);
    patches[1].computeVariableDifference(0, 0, 1, 0, max_norm_1, l2_norm_1);
    max_norm = max(max_norm_0, max_norm_1);
    l2_norm = sqrt(sqr(l2_norm_0) + sqr(l2_norm_1));

    cout << t/time << "  dt: " << dt << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm;
    cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
    ++iter;
  }

  stopTiming();
  cout << iter << " iterations" << endl;
}

#endif // TWO_PATCHES_1
