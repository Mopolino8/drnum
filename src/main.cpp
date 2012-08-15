#include "blockcfd.h"
#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/immersedboundaryreconstruction.h"
#include "fluxes/ausmplus.h"
#include "fluxes/ausmdv.h"
#include "fluxes/ausm.h"
#include "fluxes/compressiblewallflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "iterators/cartesianstandarditerator.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "iterators/cartesiandirectionalpatchoperation.h"
#include "rungekutta.h"
#include "patchgrid.h"
#include "perfectgas.h"
#include "shapes/sphere.h"
#include "shapes/halfspace.h"
#include "boundary_conditions/compressibleeulerwall.h"
#include "compressiblevariables.h"

template <typename TShape>
class MyFlux
{

  typedef ImmersedBoundaryReconstruction<Upwind2<SecondOrder>, TShape, CompressibleEulerWall> reconstruction_t;
  //typedef ImmersedBoundaryReconstruction<Upwind2<MinMod>, TShape, CompressibleEulerWall> reconstruction_t;
  //typedef ImmersedBoundaryReconstruction<Upwind1, TShape, CompressibleEulerWall> reconstruction_t;
  typedef AusmDV<reconstruction_t, PerfectGas> euler_t;
  typedef CompressibleWallFlux<reconstruction_t, PerfectGas> wall_t;
  typedef CompressibleFarfieldFlux<reconstruction_t, PerfectGas> farfield_t;

  reconstruction_t* m_Reconstruction;
  euler_t*          m_EulerFlux;
  wall_t*           m_WallFlux;
  farfield_t*       m_FarFlux;
  TShape            m_Shape;


public: // methods

  MyFlux(real u, TShape shape);
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

};


template <typename TShape>
MyFlux<TShape>::MyFlux(real u, TShape shape)
{
  m_Shape = shape;
  m_Reconstruction = new reconstruction_t(&m_Shape);
  m_EulerFlux = new euler_t(m_Reconstruction);
  m_FarFlux = new farfield_t(m_Reconstruction);
  m_FarFlux->setFarfield(1e5, 300, u, 0, 0);
  m_WallFlux = new wall_t(m_Reconstruction);
}

template <typename TShape>
inline void MyFlux<TShape>::xField
(
  CartesianPatch *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux->xField(patch, i, j, k, x, y, z, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::yField
(
  CartesianPatch *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux->yField(patch, i, j, k, x, y, z, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::zField
(
  CartesianPatch *patch,
  size_t i, size_t j, size_t k,
  real x, real y, real z,
  real A, real* flux
)
{
  m_EulerFlux->zField(patch, i, j, k, x, y, z, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarFlux->xWallP(P, i, j, k, x, y, z, A, flux);
  //m_WallFlux->xWallP(P, i, j, k, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarFlux->yWallP(P, i, j, k, x, y, z, A, flux);
  m_WallFlux->yWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarFlux->zWallP(P, i, j, k, x, y, z, A, flux);
  m_WallFlux->zWallP(P, i, j, k, x, y, z, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarFlux->xWallM(P, i, j, k, x, y, z, A, flux);
  //m_WallFlux->xWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarFlux->yWallM(P, i, j, k, x, y, z, A, flux);
  m_WallFlux->yWallM(P, i, j, k, x, y, z, A, flux);
}

template <typename TShape>
inline void MyFlux<TShape>::zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  //m_FarFlux->zWallM(P, i, j, k, x, y, z, A, flux);
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

void main1()
{
  // Testing transformations (debug only)
  CoordTransformVV ct;
  vec3_t b(1., 0., 0.);
  ct.setVector(b);
  vec3_t xyzo(0., 0., 0.);
  vec3_t xyz = ct.transform(xyzo);
  vec3_t xyzo_1 = ct.transformReverse(xyz);

  // Testing patchGrids
  PatchGrid patchGrid;
  patchGrid.setInterpolateData();
  patchGrid.setTransferPadded();

  // Define and insert some patches
  CartesianPatch patch_0;
  CartesianPatch patch_1;
  CartesianPatch patch_2;
  patch_0.setupAligned(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
  patch_0.resize(10,10,10);
  patch_1.setupAligned(8.0, 0.0, 0.0, 18.0, 10.0, 10.0);
  patch_1.resize(10,10,10);
  patch_2.setupAligned(4.0, 9.0, 0.0, 14.0, 19.0, 10.0);
  // patch_2.setupAligned(8.0, 0.0, 0.0, 18.0, 10.0, 10.0); // identical to patch_1
  patch_2.resize(10,10,10);
  patchGrid.insertPatch(&patch_0);
  patchGrid.insertPatch(&patch_1);
  patchGrid.insertPatch(&patch_2);

  // Compute dependencies
  patchGrid.computeDependencies(true);

//  patch_0.insertNeighbour(&patch_1);  /// will be automatic later
//  patch_0.insertNeighbour(&patch_2);  /// will be automatic later
//  patch_1.insertNeighbour(&patch_0);  /// will be automatic later
//  patch_1.insertNeighbour(&patch_2);  /// will be automatic later
//  patch_2.insertNeighbour(&patch_0);  /// will be automatic later
//  patch_2.insertNeighbour(&patch_1);  /// will be automatic later
//  patch_0.setTransferPadded();        /// will be envoqued by PatchGrid later
//  patch_1.setTransferPadded();        /// will be envoqued by PatchGrid later
//  patch_2.setTransferPadded();        /// will be envoqued by PatchGrid later
//  patch_0.finalizeDependencies();     /// will be automatic later
//  patch_1.finalizeDependencies();     /// will be automatic later
//  patch_2.finalizeDependencies();     /// will be automatic later
}


#define N  50
#define NI 4*N
#define NJ 1*N
#define NK 1*N

void main2()
{
  CartesianPatch patch;
  patch.setNumberOfFields(2);
  patch.setNumberOfVariables(5);
  patch.setupAligned(-3, 0, 0, 3, 1.5, 1.5);
  patch.resize(NI, NJ, NK);
  real init_var[5];
  real zero_var[5];

  real Ma = 0.3;
  real u  = Ma*sqrt(1.4*287*300.0);

  PerfectGas::primitiveToConservative(1e5, 300, u, 0, 0, init_var);
  PerfectGas::primitiveToConservative(1e5, 300, 0, 0, 0, zero_var);
  for (size_t i = 0; i < NI; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        patch.setVar(0, i, j, k, init_var);
      }
    }
  }

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.50);
  runge_kutta.addAlpha(1.00);

  Sphere shape;
  shape.setCentre(0, 0, 0);
  shape.setRadius(0.25);

  //HalfSpace shape;
  //shape.setPoint(0, 0.2, 0);
  //shape.setNormal(-1,20,0);

  shape.transform(patch.getTransformation());

  typedef Sphere shape_t;

  MyFlux<shape_t> flux(u, shape);
  CartesianDirectionalPatchOperation<5, MyFlux<shape_t> > operation(&patch, &flux);
  CartesianStandardIterator iterator(&operation);
  runge_kutta.addIterator(&iterator);

  real dt             = 3e-5;
  real dt_max         = 5e-5;
  real dt_ramp        = 1.0;
  real t_write        = 0;
  real write_interval = 1e-3;
  real total_time     = 1;

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
    real x = 0.5*patch.dx();
    for (size_t i = 0; i < NI; ++i) {
      real y = 0.5*patch.dy();
      for (size_t j = 0; j < NJ; ++j) {
        real z = 0.5*patch.dz();
        for (size_t k = 0; k < NK; ++k) {
          if (!shape.isInside(x, y, z)) {
            real p, u, v, w, T, var[5];
            patch.getVar(0, i, j, k, var);
            PerfectGas::conservativeToPrimitive(var, p, T, u, v, w);
            real a = sqrt(PerfectGas::gamma()*PerfectGas::R()*T);
            CFL_max = max(CFL_max, fabs(u)*dt/patch.dx());
            CFL_max = max(CFL_max, fabs(u+a)*dt/patch.dx());
            CFL_max = max(CFL_max, fabs(u-a)*dt/patch.dx());
            countSqrts(1);
            countFlops(10);
          } else {
            patch.setVar(0, i, j, k, zero_var);
          }
          z += patch.dz();
        }
        y += patch.dy();
      }
      x += patch.dx();
    }
    t += dt;
    t_write += dt;
    if (t_write >= write_interval) {
      ++count;
      write(patch, "testrun", count);
      t_write -= write_interval;
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
}

int main()
{
  main1();
}

