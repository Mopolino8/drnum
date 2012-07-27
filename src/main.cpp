#include "compressiblecartesianpatch.h"
#include "ausmtools.h"

#define NI 200
#define NJ 50
#define NK 50

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "fluxes/ausmplus.h"
#include "fluxes/ausmdv.h"
#include "fluxes/ausm.h"
#include "fluxes/compressibleflux.h"
#include "iterators/cartesianstandarditerator.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "rungekutta.h"

template <class TReconstruction>
class TestFlux
{

  //Ausm<TReconstruction>             m_EulerFlux;
  AusmPlus<TReconstruction>         m_EulerFlux;
  //AusmDV<TReconstruction>           m_EulerFlux;
  CompressibleFlux<TReconstruction> m_WallFlux;

public: // methods

  void x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);

  void xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);
  void zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux);

};

template <class TReconstruction>
inline void TestFlux<TReconstruction>::x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_EulerFlux.x(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_EulerFlux.y(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_EulerFlux.z(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_WallFlux.xWallP(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_WallFlux.yWallP(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_WallFlux.zWallP(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_WallFlux.xWallM(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_WallFlux.yWallM(P, i, j, k, A, flux);
}

template <class TReconstruction>
inline void TestFlux<TReconstruction>::zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, RealVec<5> &flux)
{
  m_WallFlux.zWallM(P, i, j, k, A, flux);
}

void write(CompressibleCartesianPatch &patch, QString file_name, int count)
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
  CompressibleCartesianPatch patch;
  patch.setupAligned(0, 0, 0, 8, 2, 2);
  patch.resize(NI, NJ, NK);
  for (size_t i = 0; i < NI; ++i) {
    real x = -4 + 0.5*patch.dx() + i*patch.dx();
    for (size_t j = 0; j < NJ; ++j) {
      real y = -1 + 0.5*patch.dy() + j*patch.dy();
      for (size_t k = 0; k < NK; ++k) {
        real z = -1 + 0.5*patch.dz() + k*patch.dz();
        real r = sqrt(x*x + y*y + z*z);
        if (r < 0.5) {
          patch.setState(i, j, k, 1e6, 600);
        } else {
          patch.setState(i, j, k, 1e5, 300);
        }
      }
    }
  }

  RungeKutta runge_kutta;
  //runge_kutta.addAlpha(0.25);
  //runge_kutta.addAlpha(0.50);
  runge_kutta.addAlpha(1.00);

  CartesianStandardPatchOperation<5,TestFlux<Upwind2<VanAlbada2> > > flux(&patch);
  //CartesianStandardPatchOperation<5,TestFlux<Upwind1> > flux(&patch);
  CartesianStandardIterator iterator(&flux);
  runge_kutta.addIterator(&iterator);

  real dt             = 2e-5;
  real t_write        = 0;
  real write_interval = 1e-4;
  real total_time     = 1e-3;

  int count = 0;
  int iter = 0;
  real t = 0;
  write(patch, "shock_tube", count);

  startTiming();

  while (t < total_time) {
    runge_kutta(dt);
    real CFL_max = 0;
    for (size_t i = 0; i < NI; ++i) {
      for (size_t j = 0; j < NJ; ++j) {
        for (size_t k = 0; k < NK; ++k) {
          real p, u, v, w, T;
          patch.getState(i, j, k, p, u, v, w, T);
          real a = sqrt(patch.gasGamma()*patch.gasR()*T);
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
      write(patch, "shock_tube", count);
      t_write = 0;
    }
    real max_norm, l2_norm;
    patch.computeVariableDifference(0, 0, 1, 0, max_norm, l2_norm);
    cout << t << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm << endl;
    ++iter;
  }
  write(patch, "shock_tube", -1);

  stopTiming();

  cout << iter << " iterations" << endl;
  cout << 1e-9*iter*3*15*sizeof(real)*NI*NJ*NK/(time(NULL) - global_start_time) << "GBs" << endl;

  return 0;
}

