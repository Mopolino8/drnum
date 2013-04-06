#ifndef GPUDUCT_H
#define GPUDUCT_H

#include "reconstruction/upwind1.h"
#include "reconstruction/upwind2.h"
#include "reconstruction/vanalbada.h"
#include "reconstruction/upwind2.h"
#include "fluxes/ausmplus.h"
#include "fluxes/ausmdv.h"
#include "fluxes/knp.h"
#include "fluxes/compressiblewallflux.h"
#include "fluxes/compressiblefarfieldflux.h"
#include "fluxes/compressibleviscflux.h"
#include "perfectgas.h"
#include "compressiblevariables.h"
#include "rungekutta.h"
#include "iterators/cartesianstandardpatchoperation.h"
#include "iterators/cartesianstandarditerator.h"

class MyFlux
{

  typedef Upwind2<VanAlbada> reconstruction2_t;
  typedef Upwind1            reconstruction1_t;
  typedef KNP<reconstruction1_t, PerfectGas> euler1_t;
  typedef AusmDV<reconstruction2_t, PerfectGas> euler2_t;
  typedef CompressibleSlipFlux<reconstruction1_t, PerfectGas> wall_t;
  typedef CompressibleFarfieldFlux<reconstruction1_t, PerfectGas> farfield_t;
  typedef CompressibleViscFlux<PerfectGas> viscous_t;

  reconstruction1_t*    m_Reconstruction1;
  reconstruction2_t*    m_Reconstruction2;
  euler1_t*             m_EulerFlux1;
  euler2_t*             m_EulerFlux2;
  viscous_t*            m_ViscFlux;
  wall_t*               m_WallFlux;
  farfield_t*           m_FarFlux;
  bool                  m_SecondOrder;


public: // methods

  MyFlux();

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


MyFlux::MyFlux()
{
  m_Reconstruction1 = new reconstruction1_t();
  m_Reconstruction2 = new reconstruction2_t();
  m_EulerFlux1 = new euler1_t(m_Reconstruction1);
  m_EulerFlux2 = new euler2_t(m_Reconstruction2);
  m_ViscFlux = new viscous_t();
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
}

inline void MyFlux::yWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
}

inline void MyFlux::zWallP(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
}

inline void MyFlux::xWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
}

inline void MyFlux::yWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
}

inline void MyFlux::zWallM(CartesianPatch *P, size_t i, size_t j, size_t k, real x, real y, real z, real A, real* flux)
{
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
  real delta_p;

  virtual void operator()()
  {
    real p, T;
    vec3_t U;
    real var[5];
    for (size_t i = 0; i < patch->sizeI(); ++i) {
      for (size_t j = 0; j < patch->sizeJ(); ++j) {

        // bottom
        patch->getVar(0, i, j, 1, var);
        var[1] = -var[1];
        var[2] = -var[2];
        var[3] = -var[3];
        patch->setVar(0, i, j, 0, var);

        // top
        patch->getVar(0, i, j, patch->sizeK() - 2, var);
        var[1] = -var[1];
        var[2] = -var[2];
        var[3] = -var[3];
        patch->setVar(0, i, j, patch->sizeK() - 1, var);

      }
    }
    for (size_t j = 0; j < patch->sizeJ(); ++j) {
      for (size_t k = 0; k < patch->sizeK(); ++k) {

        // left
        patch->getVar(0, patch->sizeI() - 2, j, k, var);
        PerfectGas::conservativeToPrimitive(var, p, T, U);
        p = p + 50;
        PerfectGas::primitiveToConservative(p, T, U, var);
        patch->setVar(0, 0, j, k, var);

        // right
        patch->getVar(0, 1, j, k, var);
        PerfectGas::conservativeToPrimitive(var, p, T, U);
        p = p - 50;
        PerfectGas::primitiveToConservative(p, T, U, var);
        patch->setVar(0, patch->sizeI() - 1, j, k, var);

      }
    }
    for (size_t i = 0; i < patch->sizeI(); ++i) {
      for (size_t k = 0; k < patch->sizeK(); ++k) {

        // front
        patch->getVar(0, i, patch->sizeJ() - 2, k, var);
        patch->setVar(0, i, 0, k, var);

        // back
        patch->getVar(0, i, 1, k, var);
        patch->setVar(0, i, patch->sizeJ() - 1, k, var);

      }
    }
  }
};

void run()
{
  CartesianPatch patch;
  patch.setNumberOfFields(2);
  patch.setNumberOfVariables(5);
  real height = 0.01;
  real length = 2*height;
  patch.setupAligned(0, 0 , 0, length, 0.001, height);
  size_t N = 10;
  size_t NI = 20*N;
  size_t NJ = N;
  size_t NK = 10*N;
  patch.resize(NI, NJ, NK);
  real init_var[5];

  real Ma = 0.3;
  real p  = 1e5;
  real Dp = 106;
  real T  = 300;
  real u  = Ma*sqrt(PerfectGas::gamma()*PerfectGas::R()*T);

  for (size_t i = 0; i < NI; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {
        real x = 0.5*patch.dx() + i*patch.dx();
        real var[5];
        PerfectGas::primitiveToConservative(p - x*Dp/length, T, vec3_t(u,0,0), var);
        patch.setVar(0, i, j, k, var);
      }
    }
  }

  BC bc;
  bc.patch = &patch;
  bc.delta_p = Dp;

  RungeKutta runge_kutta;
  runge_kutta.addPostOperation(&bc);
  runge_kutta.addAlpha(0.25);
  runge_kutta.addAlpha(0.5);
  runge_kutta.addAlpha(1.000);

  MyFlux flux;
  CartesianStandardPatchOperation<5, MyFlux> operation(&patch, &flux);
  CartesianStandardIterator iterator(&operation);
  runge_kutta.addIterator(&iterator);

  real time           = NI*patch.dx()/u;
  real dt             = 2e-3*patch.dx();
  real dt2            = dt;
  real t_switch       = 0*time;
  real t_write        = 0;
  real write_interval = 0.1*time;
  real total_time     = 4*time;

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
    cout << t/time << "  dt: " << dt << "  CFL: " << CFL_max << "  max: " << max_norm << "  L2: " << l2_norm;
    cout << "  min(rho): " << rho_min << "  max(rho): " << rho_max << endl;
    ++iter;
  }

  stopTiming();
  cout << iter << " iterations" << endl;
}

#endif // GPUDUCT_H
