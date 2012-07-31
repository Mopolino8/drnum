#define NI 50
#define NJ 50
#define NK 50

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
#include "perfectgas.h"
#include "shapes/sphere.h"
#include "boundary_conditions/compressibleeulerwall.h"


typedef ImmersedBoundaryReconstruction<Upwind2<VanAlbada2>, Sphere, CompressibleEulerWall> reconstruction_t;
//typedef Upwind2<VanAlbada2> reconstruction_t;
typedef AusmDV<reconstruction_t, PerfectGas> euler_t;
typedef CompressibleWallFlux<reconstruction_t, PerfectGas> wall_t;
typedef CompressibleFarfieldFlux<reconstruction_t, PerfectGas> farfield_t;

class MyFlux
{

  reconstruction_t* m_Reconstruction;
  euler_t*          m_EulerFlux;
  wall_t*           m_WallFlux;
  farfield_t*       m_FarFlux;

  Sphere            m_Sphere;


public: // methods

  MyFlux(real u);

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


MyFlux::MyFlux(real u)
{
  m_Sphere.setCentre(0, 0, 0);
  m_Sphere.setRadius(1.0);
  m_Reconstruction = new reconstruction_t(&m_Sphere);
  //m_Reconstruction = new reconstruction_t();
  m_EulerFlux = new euler_t(m_Reconstruction);
  m_FarFlux = new farfield_t(m_Reconstruction);
  m_FarFlux->setFarfield(1e5, 300, u, 0, 0);
  m_WallFlux = new wall_t(m_Reconstruction);
}

inline void MyFlux::x(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_EulerFlux->x(P, i, j, k, A, flux);
}

inline void MyFlux::y(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_EulerFlux->y(P, i, j, k, A, flux);
}

inline void MyFlux::z(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_EulerFlux->z(P, i, j, k, A, flux);
}

inline void MyFlux::xWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux->xWallP(P, i, j, k, A, flux);
  //m_WallFlux->xWallP(P, i, j, k, A, flux);
}

inline void MyFlux::yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  //m_FarFlux->yWallP(P, i, j, k, A, flux);
  m_WallFlux->yWallP(P, i, j, k, A, flux);
}

inline void MyFlux::zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  //m_FarFlux->zWallP(P, i, j, k, A, flux);
  m_WallFlux->zWallP(P, i, j, k, A, flux);
}

inline void MyFlux::xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  m_FarFlux->xWallM(P, i, j, k, A, flux);
  //m_WallFlux->xWallM(P, i, j, k, A, flux);
}

inline void MyFlux::yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  //m_FarFlux->yWallM(P, i, j, k, A, flux);
  m_WallFlux->yWallM(P, i, j, k, A, flux);
}

inline void MyFlux::zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real A, real* flux)
{
  //m_FarFlux->zWallM(P, i, j, k, A, flux);
  m_WallFlux->zWallM(P, i, j, k, A, flux);
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
  patch.setupAligned(-5, 0, 0, 5, 5, 5);
  patch.resize(NI, NJ, NK);
  for (size_t i = 0; i < NI; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        real var[5];
        PerfectGas::primitiveToConservative(1e5, 300, 100, 0, 0, var);
        patch.setVar(0, i, j, k, var);
      }
    }
  }

  RungeKutta runge_kutta;
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.50);
  runge_kutta.addAlpha(1.00);

  MyFlux flux(100.0);
  CartesianDirectionalPatchOperation<5, MyFlux> operation(&patch, &flux);
  CartesianStandardIterator iterator(&operation);
  runge_kutta.addIterator(&iterator);

  real dt             = 1e-4;
  real dt_max         = 1e-4;
  real dt_ramp        = 1.1;
  real t_write        = 0;
  real write_interval = 1e-3;
  real total_time     = 1e-3;

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

