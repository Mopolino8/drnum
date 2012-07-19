#include "compressiblecartesianpatch.h"
#include "ausmtools.h"

#define NI 1000
#define NJ 10
#define NK 10

#include "reconstruction/first_order.h"

class TestPatch : public CompressibleCartesianPatch
{

public:

  virtual void preStep();
  virtual void postStep();
  void correctBoundaries();

};

inline void TestPatch::preStep()
{
  COMPR_NEW_VARIABLES;
  COMPR_RES_VARIABLES;

  for (size_t i = 1; i < NI - 3; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {

        #include "fluxes/ausm_plus_x.h"

      }
    }
  }
}

inline void TestPatch::postStep()
{
  for (size_t j = 0; j < NJ; ++j) {
    for (size_t k = 0; k < NK; ++k) {
      for (size_t i_var = 0; i_var < 5; ++i_var) {
        f(getVariable(i_new, i_var), 0, j, k) = f(getVariable(i_new, i_var), 2, j, k);
        f(getVariable(i_new, i_var), 1, j, k) = f(getVariable(i_new, i_var), 2, j, k);
        f(getVariable(i_new, i_var), NI-2, j, k) = f(getVariable(i_new, i_var), NI-3, j, k);
        f(getVariable(i_new, i_var), NI-1, j, k) = f(getVariable(i_new, i_var), NI-3, j, k);
      }
    }
  }
}

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
  real dt = 1e-5;
  real t = 0;
  real alpha[3] = {0.25, 0.5, 1};
  int count = 0;
  write(P, count);
  while (t < 1e-3) {
    P.copyField(P.i_new, P.i_old);
    for (int i_rk = 0; i_rk < 3; ++i_rk) {
      P.subStep(dt*alpha[i_rk]);
    }
    t += dt;
    ++count;
    write(P, count);
    cout << t << endl;
  }
}

