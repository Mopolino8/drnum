#include "compressiblecartesianpatch.h"
#include "ausmtools.h"

#define NI 1000
#define NJ 5
#define NK 5

#include "reconstruction/upwind2.h"
#include "fluxes/ausmplus.h"
#include "iterators/cartesianstandarditerator.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "rungekutta.h"
#include "patchgrid.h"

template <class TReconstruction>
class TestFlux
{

  AusmPlus<TReconstruction>         m_EulerFlux;
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


/*
typedef CompressibleCartesianPatch<AusmPlus<Upwind2<VanAlbada2> > > TestPatch;
//typedef CompressibleCartesianPatch<AusmPlus<Upwind2<FirstOrder> > > TestPatch;

void write(TestPatch &P, int count)
{
  QString file_name;
  file_name.setNum(count);
  while (file_name.size() < 6) {
    file_name = "0" + file_name;
  }
  file_name = "shock_tube_" + file_name;
  P.writeToVtk(file_name);
}

int main()
{
  TestPatch P;
  P.setupAligned(0, 0, 0, 10.0, 0.1, 0.1);
  P.resize(NI, NJ, NK);
  for (size_t i = 0; i < NI; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        if (i < NI/2) {
          P.setState(i, j, k, 1e6, 300);
        } else {
          P.setState(i, j, k, 1e5, 300);
        }
      }
    }
  }
  real t = 0;
  real alpha[3] = {0.25, 0.5, 1};
  int count = 0;

  real dt             = 1e-5;
  real t_write        = 0;
  real write_interval = 1e-4;
  real total_time     = 5e-3;

  write(P, count);
  while (t < total_time) {
    P.copyField(P.i_new, P.i_old);
    for (int i_rk = 0; i_rk < 3; ++i_rk) {
      P.subStep(dt*alpha[i_rk]);
    }
    real CFL_max = 0;
    for (size_t i = 0; i < NI; ++i) {
      for (size_t j = 0; j < NJ; ++j) {
        for (size_t k = 0; k < NK; ++k) {
          real p, u, v, w, T;
          P.getState(i, j, k, p, u, v, w, T);
          real a = sqrt(P.gasGamma()*P.gasR()*T);
          CFL_max = max(CFL_max, fabs(u)*dt/P.dx());
          CFL_max = max(CFL_max, fabs(u+a)*dt/P.dx());
          CFL_max = max(CFL_max, fabs(u-a)*dt/P.dx());
        }
      }
    }
    t += dt;
    t_write += dt;
    if (t_write >= write_interval) {
      ++count;
      write(P, count);
      t_write = 0;
    }
    real max_norm, l2_norm;
    P.computeVariableDifference(P.i_new, 0, P.i_old, 0, max_norm, l2_norm);
    cout << t << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm << endl;
  }
}
*/

void write(CompressibleCartesianPatch &patch, int count)
{
  QString file_name;
  file_name.setNum(count);
  while (file_name.size() < 6) {
    file_name = "0" + file_name;
  }
  file_name = "shock_tube_" + file_name;
  patch.writeToVtk(file_name);
}

int main()
{
  /// Testing transformations (debug only)
  CoordTransformVV ct;
  vec3_t b(1., 0., 0.);
  ct.setVector(b);
  vec3_t xyzo(0., 0., 0.);
  vec3_t xyz = ct.transform(xyzo);
  vec3_t xyzo_1 = ct.transformReverse(xyz);

  /// Testing patchGrids
  PatchGrid patchGrid;
  patchGrid.initLists(10,10);
  CartesianPatch patch_0;
  CartesianPatch patch_1;
  patch_0.setupAligned( 0., 0., 0., 10.0, 10.0, 10.0);
  patch_0.resize(10,10,10);
  patch_0.setInterpolateData();
  patch_1.setupAligned(10., 0., 0., 20.0, 10.0, 10.0);
  patch_1.resize(10,10,10);
  patch_1.setInterpolateData();
  patchGrid.insertPatch(&patch_0);
  patchGrid.insertPatch(&patch_1);
  // patchGrid.computeDependencies(true);
  patch_0.insertNeighbour(&patch_1);  /// will be automatic later
  patch_1.insertNeighbour(&patch_0);  /// will be automatic later




  CompressibleCartesianPatch patch; //@@ kann ersetzt werden
  patch.setupAligned(0, 0, 0, 10.0, 0.1, 0.1);
  patch.resize(NI, NJ, NK);
  for (size_t i = 0; i < NI; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        if (i < NI/2) {
          patch.setState(i, j, k, 1e6, 300);
        } else {
          patch.setState(i, j, k, 1e5, 300);
        }
      }
    }
  }

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.50);
  runge_kutta.addAlpha(1.00);

  CartesianStandardPatchOperation<5,TestFlux<Upwind2<VanAlbada2> > > flux(&patch);
  CartesianStandardIterator iterator(&flux);
  runge_kutta.addIterator(&iterator);

  real dt             = 1e-5;
  real t_write        = 0;
  real write_interval = 1e-4;
  real total_time     = 5e-3;

  int count = 0;
  real t = 0;
  write(patch, count);
  while (t < total_time) {
    runge_kutta(dt);

    real CFL_max = 0;
    /*
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
      write(patch, count);
      t_write = 0;
    }
    */
    real max_norm, l2_norm;
    patch.computeVariableDifference(0, 0, 1, 0, max_norm, l2_norm);
    cout << t << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm << endl;
  }

  return 0;
}

