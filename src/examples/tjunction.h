#ifndef TJUNCTION_H
#define TJUNCTION_H

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/minmod.h"
#include "reconstruction/secondorder.h"
#include "fluxes/knp.h"
#include "fluxes/vanleer.h"
#include "fluxes/ausmplus.h"
#include "fluxes/compressiblewallflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "boundary_conditions/compressibleeulerwall.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "iterators/cartesianstandarditerator.h"

class MyFlux
{

  typedef Upwind1 reconstruction1_t;
  typedef Upwind2<VanAlbada> reconstruction2_t;
  typedef VanLeer<reconstruction1_t, PerfectGas> euler1_t;
  typedef AusmPlus<reconstruction2_t, PerfectGas> euler2_t;
  typedef CompressibleWallFlux<reconstruction1_t, PerfectGas> wall_t;
  typedef CompressibleFarfieldFlux<reconstruction1_t, PerfectGas> farfield_t;
  typedef CompressibleViscFlux<PerfectGas> viscous_t;

  reconstruction1_t* m_Reconstruction1;
  reconstruction2_t* m_Reconstruction2;
  euler1_t*          m_EulerFlux1;
  euler2_t*          m_EulerFlux2;
  viscous_t*         m_ViscFlux;
  wall_t*            m_WallFlux;
  farfield_t*        m_FarFlux;
  bool               m_SecondOrder;


public: // methods

  MyFlux(real p, real T, real u);
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


MyFlux::MyFlux(real p, real T, real u)
{
  m_Reconstruction1 = new reconstruction1_t();
  m_Reconstruction2 = new reconstruction2_t();
  m_EulerFlux1 = new euler1_t(m_Reconstruction1);
  m_EulerFlux2 = new euler2_t(m_Reconstruction2);
  m_ViscFlux = new viscous_t();
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
  m_ViscFlux->xField(patch, i, j, k, x, y, z, A, flux);
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
  m_ViscFlux->yField(patch, i, j, k, x, y, z, A, flux);
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
  m_ViscFlux->zField(patch, i, j, k, x, y, z, A, flux);
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


struct BC : public GenericOperation
{
  CartesianPatch *patch;
  real* init_var1;
  real* init_var2;
  size_t N;

  virtual void operator()()
  {
    real var[5];
    real p, T, p1, T1, p2, T2;
    vec3_t U, U1, U2;
    PerfectGas::conservativeToPrimitive(init_var1, p1, T1, U1);
    PerfectGas::conservativeToPrimitive(init_var2, p2, T2, U2);

    // lower wall
    for (size_t i = 0; i < patch->sizeI(); ++i) {
      patch->getVar(0, i, 0, 1, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      U[2] = -U[2];
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, i, 0, 0, var);
    }

    // upper wall (left)
    for (size_t i = 0; i < 2*N; ++i) {
      patch->getVar(0, i, 0, N-1, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      U[2] = -U[2];
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, i, 0, N, var);
    }

    // upper wall (right)
    for (size_t i = 3*N+1; i < patch->sizeI(); ++i) {
      patch->getVar(0, i, 0, N-1, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      U[2] = -U[2];
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, i, 0, N, var);
    }

    // left wall
    for (size_t k = N+1; k < patch->sizeK(); ++k) {
      patch->getVar(0, 2*N+1, 0, k, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      U[0] = -U[0];
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, 2*N, 0, k, var);
    }

    // right wall
    for (size_t k = N+1; k < patch->sizeK(); ++k) {
      patch->getVar(0, 3*N-1, 0, k, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      U[0] = -U[0];
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, 3*N, 0, k, var);
    }

    // left corner
    patch->getVar(0, 2*N+1, 0, N-1, var);
    PerfectGas::conservativeToPrimitive(var, p, T, U);
    PerfectGas::primitiveToConservative(p, T, var);
    patch->setVar(0, 2*N, 0, N, var);

    // right corner
    patch->getVar(0, 3*N-1, 0, N-1, var);
    PerfectGas::conservativeToPrimitive(var, p, T, U);
    PerfectGas::primitiveToConservative(p, T, var);
    patch->setVar(0, 3*N, 0, N, var);

    // left inlet
    for (size_t k = 0; k <= N; ++k) {
      patch->getVar(0, 1, 0, k, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      T = T1;
      U = U1;
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, 0, 0, k, var);
    }

    // right inlet
    for (size_t k = 0; k <= N; ++k) {
      patch->getVar(0, patch->sizeI()-2, 0, k, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      T = T2;
      U = U2;
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, patch->sizeI()-1, 0, k, var);
    }

    // outlet
    for (size_t i = 2*N+1; i < 3*N; ++i) {
      patch->getVar(0, i, 0, patch->sizeK()-2, var);
      PerfectGas::conservativeToPrimitive(var, p, T, U);
      p = p1;
      PerfectGas::primitiveToConservative(p, T, U, var);
      patch->setVar(0, i, 0, patch->sizeK()-1, var);
    }

    // dead zones
    for (size_t i = 0; i < 2*N; ++i) {
      for (size_t k = N+1; k < patch->sizeK(); ++k) {
        PerfectGas::primitiveToConservative(p1, 0.5*(T1+T2), var);
        patch->setVar(0, i, 0, k, var);
      }
    }
    for (size_t i = 3*N+1; i < patch->sizeI(); ++i) {
      for (size_t k = N+1; k < patch->sizeK(); ++k) {
        PerfectGas::primitiveToConservative(p1, 0.5*(T1+T2), var);
        patch->setVar(0, i, 0, k, var);
      }
    }

  }

};


void run()
{
  size_t N = 50;
  real Ma1 = 0.05;
  real p1  = 1e5;
  real T1  = 600;
  real u1  = Ma1*sqrt(PerfectGas::gamma()*PerfectGas::R()*T1);
  //real Ma2 = 0.1;
  real p2  = 1e5;
  real T2  = 300;
  real u2  = -u1;//Ma2*sqrt(PerfectGas::gamma()*PerfectGas::R()*T2);
  real Re  = 1000;
  real CFL = 0.7;
  real L   = PerfectGas::mu()*Re/u1;

  cout << "L = " << L << endl;

  CartesianPatch patch;
  patch.setNumberOfFields(2);
  patch.setNumberOfVariables(5);
  patch.setupAligned(-2.5*L, -0.01*L, 0, 2.5*L, 0.01*L, 5*L);
  size_t NI = 5*N;
  size_t NJ = 1;
  size_t NK = 5*N;
  patch.resize(NI, NJ, NK);
  real init0_var[5];
  real init1_var[5];
  real init2_var[5];

  PerfectGas::primitiveToConservative(p2, T2, 0, 0, 0, init0_var);
  PerfectGas::primitiveToConservative(p1, T1, u1, 0, 0, init1_var);
  PerfectGas::primitiveToConservative(p2, T2, u2, 0, 0, init2_var);
  for (size_t i = 0; i < NI; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        patch.setVar(0, i, j, k, init0_var);
      }
    }
  }

  BC bc;
  bc.patch = &patch;
  bc.init_var1 = init1_var;
  bc.init_var2 = init2_var;
  bc.N = N;

  RungeKutta runge_kutta;
  runge_kutta.addPostOperation(&bc);
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);


  MyFlux flux(p1, T1, u1);
  CartesianStandardPatchOperation<5, MyFlux> operation(&patch, &flux);
  CartesianStandardIterator iterator(&operation);
  iterator.setI1(1);
  runge_kutta.addIterator(&iterator);

  real time           = L/(u1 + sqrt(PerfectGas::gamma()*PerfectGas::R()*T1));
  real dt             = CFL*time/N;
  real dt2            = dt;
  real t_switch       = 0*L/u1;
  real t_write        = 0;
  real write_interval = L/u1/10;
  real total_time     = 100*L/u1;

  int count = 0;
  int iter = 0;
  real t = 0;
  write(patch, "testrun", count);

  cout << "Press <ENTER> to start!";
  cin.get();

  startTiming();

  while (t < total_time) {

    if (t > t_switch) {
      flux.secondOrderOn();
      dt = dt2;
    }
    runge_kutta(dt);

    real CFL_max = 0;
    real x = 0.5*patch.dx();
    real rho_min = 1000;
    real rho_max = 0;
    for (size_t i = 0; i < NI; ++i) {
      real y = 0.5*patch.dy();
      for (size_t j = 0; j < NJ; ++j) {
        real z = 0.5*patch.dz();
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
    cout << u1*t/L << "  dt: " << dt << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm;
    cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
    ++iter;
  }

  stopTiming();
  cout << iter << " iterations" << endl;
}

#endif // TJUNCTION_H
