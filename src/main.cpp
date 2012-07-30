#define NI 200
#define NJ 50
#define NK 50

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "fluxes/ausmplus.h"
#include "fluxes/ausmdv.h"
#include "fluxes/ausm.h"
#include "fluxes/compressiblewallflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "iterators/cartesianstandarditerator.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "iterators/cartesiandirectionalpatchoperation.h"
#include "rungekutta.h"
#include "perfectgas.h"

template <typename TReconstruction>
class TestFlux
{

  AusmDV<TReconstruction, PerfectGas> m_EulerFlux;
  CompressibleWallFlux<TReconstruction, PerfectGas> m_WallFlux;
  CompressibleFarfieldFlux<TReconstruction, PerfectGas> m_FarFlux;

public: // methods

  TestFlux(real u);

  void x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);
  void y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);
  void z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);

  void xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);
  void yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);
  void zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);
  void xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);
  void yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);
  void zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux);

};


template <typename TReconstruction>
TestFlux<TReconstruction>::TestFlux(real u)
{
  m_FarFlux.setFarfield(1e5, 300, u, 0, 0);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_EulerFlux.x(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_EulerFlux.y(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_EulerFlux.z(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux.xWallP(P, i, j, k, A, flux);
  //m_WallFlux.xWallP(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux.yWallP(P, i, j, k, A, flux);
  //m_WallFlux.yWallP(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux.zWallP(P, i, j, k, A, flux);
  //m_WallFlux.zWallP(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux.xWallM(P, i, j, k, A, flux);
  //m_WallFlux.xWallM(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux.yWallM(P, i, j, k, A, flux);
  //m_WallFlux.yWallM(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux.zWallM(P, i, j, k, A, flux);
  //m_WallFlux.zWallM(P, i, j, k, A, flux);
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
  cout << "writing to file \"" << qPrintable(file_name) << endl;
  patch.writeToVtk(file_name);
}

int main()
{
  CartesianPatch patch;
  patch.setNumberOfFields(2);
  patch.setNumberOfVariables(5);
  patch.setupAligned(0, 0, 0, 8, 2, 2);
  patch.resize(NI, NJ, NK);
  for (size_t i = 0; i < NI; ++i) {
    real x = 0.5*patch.dx() + i*patch.dx();
    for (size_t j = 0; j < NJ; ++j) {
      real y = 0.5*patch.dy() + j*patch.dy();
      for (size_t k = 0; k < NK; ++k) {
        real z = 0.5*patch.dz() + k*patch.dz();
        real r = sqrt(x*x + y*y + z*z);
        real var[5];
        if (r < 0.5) {
          PerfectGas::primitiveToConservative(1e5, 300, 0, 0, 0, var);
        } else {
          PerfectGas::primitiveToConservative(1e5, 300, 0, 0, 0, var);
        }
        patch.setVar(0, i, j, k, var);
      }
    }
  }

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.50);
  runge_kutta.addAlpha(1.00);

  //CartesianStandardPatchOperation<5,TestFlux<Upwind2<VanAlbada2> > > flux(&patch);
  TestFlux<Upwind2<VanAlbada2> > flux(100.0);
  CartesianDirectionalPatchOperation<5,TestFlux<Upwind2<VanAlbada2> > > operation(&patch, &flux);
  //CartesianStandardPatchOperation<5,TestFlux<Upwind1> > flux(&patch);
  CartesianStandardIterator iterator(&operation);
  runge_kutta.addIterator(&iterator);

  real dt             = 1e-4;
  real dt_max         = 1e-4;
  real dt_ramp        = 1.1;
  real t_write        = 0;
  real write_interval = 2e-3;
  real total_time     = 1e-1;

  int count = 0;
  int iter = 0;
  real t = 0;
  write(patch, "testrun", count);

  cout << "Press <ENTER> to start!";
  cin.get();

  startTiming();

  while (t < total_time) {
    runge_kutta(dt);
    real CFL_max = 0;
    for (size_t i = 0; i < NI; ++i) {
      for (size_t j = 0; j < NJ; ++j) {
        for (size_t k = 0; k < NK; ++k) {
          real p, u, v, w, T, var[5];
          patch.getVar(0, i, j, k, var);
          PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
          real a = sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
          CFL_max = max(CFL_max, fabs(u)*dt/patch.dx());
          CFL_max = max(CFL_max, fabs(u+a)*dt/patch.dx());
          CFL_max = max(CFL_max, fabs(u-a)*dt/patch.dx());
        }
      }
    }
    t += dt;
    t_write += dt;
    if (t_write >= write_interval) {
      ++count;
      write(patch, "testrun", count);
      t_write = 0;
    }
    real max_norm, l2_norm;
    patch.computeVariableDifference(0, 0, 1, 0, max_norm, l2_norm);
    cout << t << "  dt: " << dt << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm << endl;
    ++iter;
    dt = min(dt_max, dt_ramp*dt);
  }

  stopTiming();
  write(patch, "testrun", -1);

  cout << iter << " iterations" << endl;
  cout << 1e-9*iter*3*15*sizeof(real)*NI*NJ*NK/(time(NULL) - global_start_time) << "GBs" << endl;

  return 0;
}

