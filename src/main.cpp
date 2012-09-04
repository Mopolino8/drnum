
#include "blockcfd.h"

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/upwindcentral.h"
#include "reconstruction/immersedboundaryreconstruction.h"
#include "reconstruction/minmod.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/vanleerlim.h"
#include "reconstruction/roelim.h"
#include "reconstruction/secondorder.h"

#include "fluxes/ausmplus.h"
#include "fluxes/ausmdv.h"
#include "fluxes/ausm.h"
#include "fluxes/kt.h"
#include "fluxes/knp.h"
#include "fluxes/roe.h"
#include "fluxes/vanleer.h"
#include "fluxes/compressiblewallflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "iterators/cartesianstandarditerator.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "iterators/cartesiandirectionalpatchoperation.h"
#include "rungekutta.h"
#include "patchgrid.h"
#include "perfectgas.h"
#include "shapes/sphere.h"
#include "shapes/halfspace.h"
#include "shapes/noshape.h"
#include "shapes/triangulatedshape.h"
#include "shapes/box.h"
#include "shapes/cylindery.h"
#include "boundary_conditions/compressibleeulerwall.h"
#include "compressiblevariables.h"

#define PROJ    Upwind2
#define FLUX    KT
#define LIMITER MinMod
#define EULER

template <typename TShape>
class MyFlux
{

  typedef ImmersedBoundaryReconstruction<PROJ<LIMITER>, TShape, CompressibleEulerWall> reconstruction2_t;
  typedef ImmersedBoundaryReconstruction<PROJ<MinMod>, TShape, CompressibleEulerWall> reconstruction1_t;
  //typedef ImmersedBoundaryReconstruction<Upwind1, TShape, CompressibleEulerWall> reconstruction1_t;
  typedef FLUX<reconstruction1_t, PerfectGas> euler1_t;
  typedef FLUX<reconstruction2_t, PerfectGas> euler2_t;
  typedef CompressibleWallFlux<reconstruction1_t, PerfectGas> wall_t;
  typedef CompressibleFarfieldFlux<reconstruction1_t, PerfectGas> farfield_t;

  reconstruction1_t*    m_Reconstruction1;
  reconstruction2_t*    m_Reconstruction2;
  euler1_t*             m_EulerFlux1;
  euler2_t*             m_EulerFlux2;
  wall_t*               m_WallFlux;
  farfield_t*           m_FarFlux;
  TShape*               m_Shape;
  CompressibleViscFlux<PerfectGas>* m_ViscousFlux;

  bool                  m_SecondOrder;


public: // methods

  MyFlux(real u, TShape* shape);
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


template <typename TShape>
MyFlux<TShape>::MyFlux(real u, TShape* shape)
{
  m_Shape = shape;
  m_Reconstruction1 = new reconstruction1_t(m_Shape);
  m_Reconstruction2 = new reconstruction2_t(m_Shape);
  m_EulerFlux1 = new euler1_t(m_Reconstruction1);
  m_EulerFlux2 = new euler2_t(m_Reconstruction2);
  m_ViscousFlux = new CompressibleViscFlux<PerfectGas>();
  m_FarFlux = new farfield_t(m_Reconstruction1);
  m_FarFlux->setFarfield(1e5, 300, u, 0, 0);
  m_WallFlux = new wall_t(m_Reconstruction1);
  m_SecondOrder = false;
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
  if (m_SecondOrder) m_EulerFlux2->xField(patch, i, j, k, x, y, z, A, flux);
  else               m_EulerFlux1->xField(patch, i, j, k, x, y, z, A, flux);
#ifndef EULER
  m_ViscousFlux->xField(patch, i, j, k, x, y, z, A, flux);
#endif
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
  if (m_SecondOrder) m_EulerFlux2->yField(patch, i, j, k, x, y, z, A, flux);
  else               m_EulerFlux1->yField(patch, i, j, k, x, y, z, A, flux);
#ifndef EULER
  m_ViscousFlux->yField(patch, i, j, k, x, y, z, A, flux);
#endif
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
  if (m_SecondOrder) m_EulerFlux2->zField(patch, i, j, k, x, y, z, A, flux);
  else               m_EulerFlux1->zField(patch, i, j, k, x, y, z, A, flux);
#ifndef EULER
  m_ViscousFlux->zField(patch, i, j, k, x, y, z, A, flux);
#endif
}

template <typename TShape>
inline void MyFlux<TShape>::xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
  m_FarFlux->xWallP(P, i, j, k, x, y, z, A, flux);
  //m_WallFlux->xWallP(P, i, j, k, x, y, z, A, flux);
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

//  // Define and insert some patches
//  CartesianPatch patch_0;
//  CartesianPatch patch_1;
//  CartesianPatch patch_2;
//  patch_0.setupAligned(0.0, 0.0, 0.0, 10.0, 10.0, 10.0);
//  patch_0.resize(10,10,10);
//  patch_1.setupAligned(8.0, 0.0, 0.0, 18.0, 10.0, 10.0);
//  patch_1.resize(10,10,10);
//  patch_2.setupAligned(4.0, 9.0, 0.0, 14.0, 19.0, 10.0);
//  // patch_2.setupAligned(8.0, 0.0, 0.0, 18.0, 10.0, 10.0); // identical to patch_1
//  patch_2.resize(10,10,10);
//  patchGrid.insertPatch(&patch_0);
//  patchGrid.insertPatch(&patch_1);
//  patchGrid.insertPatch(&patch_2);

  // Define and insert lots of patches
  size_t nn_i = 20;
  size_t nn_j = 20;
  size_t nn_k = 20;
  CartesianPatch* patches;
  patches = new CartesianPatch[nn_i*nn_j*nn_k];
  for (size_t i = 0; i < nn_i; i++) {
    for (size_t j = 0; j < nn_j; j++) {
      for (size_t k = 0; k < nn_k; k++) {
        size_t l = i*nn_j*nn_k + j*nn_k + k;
        real x_low = 5*i;
        real y_low = 5*j;
        real z_low = 5*k;
        patches[l].setupAligned(x_low, y_low, z_low, x_low+10.0, y_low+10.0, z_low+10.0);
        patches[l].resize(10,10,10);
        patchGrid.insertPatch(&patches[l]);
      }
    }
  }

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


//#define N  100
//#define NI 3*N
//#define NJ 3
//#define NK N

void entropyCorrection(CartesianPatch* patch, size_t i, size_t k)
{
  for (size_t j = 0; j < patch->sizeJ(); ++j) {
    real var0[5];
    patch->getVar(0, i, j, k, var0);
    real p0, T0;
    vec3_t U0;
    PerfectGas::conservativeToPrimitive(var0, p0, T0, U0);
    real ht = PerfectGas::cp()*T0 + 0.5*U0.abs2();
    static const int N = 6;
    size_t ic[N] = { i+1, i+2, i+3, i+4, i+1, i+2 };
    size_t kc[N] = { k,   k,   k,   k,   k+1, k+1 };
    for (int i_cell = 0; i_cell < N; ++i_cell) {
      real var[5];
      patch->getVar(0, ic[i_cell], j, kc[i_cell], var);
      real p, T;
      vec3_t U;
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      if (U.abs2() > 1) {
        T = T0*pow(p/p0, (PerfectGas::gamma() - 1)/PerfectGas::gamma());
        U[2] = 0;
        U.normalise();
        U *= sqrt(2*(ht - PerfectGas::cp()*T));
        PerfectGas::primitiveToConservative(p, T, U, var);
      }
      patch->setVar(0, ic[i_cell], j, kc[i_cell], var);
    }
  }
}

void main2()
{
  CartesianPatch patch;
  patch.setNumberOfFields(2);
  patch.setNumberOfVariables(5);
  patch.setupAligned(0, -0.01, 0, 3, 0.01, 1);
  size_t N = 100;
  size_t NI = 3*N;
  size_t NJ = 1;
  size_t NK = 1*N;
  patch.resize(NI, NJ, NK);
  real init_var[5];
  real zero_var[5];

  real Ma = 3;
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
  //runge_kutta.addAlpha(0.059);
  //runge_kutta.addAlpha(0.14);
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

  typedef Box shape_t;
  shape_t shape;
  shape.setGeometry(0.6, -1, -1, 10, 1, 0.2);

  shape.transform(patch.getTransformation());

  MyFlux<shape_t> flux(u, &shape);
  CartesianStandardPatchOperation<5, MyFlux<shape_t> > operation(&patch, &flux);
  CartesianStandardIterator iterator(&operation);
  runge_kutta.addIterator(&iterator);

  real T              = 1.4401e-3;
  real dt             = 5e-4/N;
  real dt2            = 1*dt;
  real t_switch       = 0*T;//2e-2;
  real t_write        = 0;
  real write_interval = T;
  real total_time     = 10*write_interval;

  int count = 0;
  int iter = 0;
  real t = 0;
  write(patch, "testrun", count);

  cout << "Press <ENTER> to start!";
  cin.get();

  startTiming();

  // compute "corner" indices
  size_t i_corner, k_corner;
  {
    real x = 0.6 - 0.5*patch.dx();
    real z = 0.2 + 0.5*patch.dz();
    i_corner = x/patch.dx();
    k_corner = z/patch.dz();
  }

  while (t < total_time) {
    if (t > t_switch) {
      flux.secondOrderOn();
      dt = dt2;
    }
    runge_kutta(dt);
    //entropyCorrection(&patch, i_corner, k_corner);

    real CFL_max = 0;
    real x = 0.5*patch.dx();
    real rho_min = 1000;
    real rho_max = 0;
    for (size_t i = 0; i < NI; ++i) {
      real y = 0.5*patch.dy();
      for (size_t j = 0; j < NJ; ++j) {
        real z = 0.5*patch.dz();
        for (size_t k = 0; k < NK; ++k) {
          if (!shape.isInside(x, y, z)) {
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
    cout << t/T << "  dt: " << dt << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm;
    cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
    ++iter;
  }

  stopTiming();
  write(patch, "testrun", -1);

  cout << iter << " iterations" << endl;
}

int main()
{
  //main1();
  main2();
}

