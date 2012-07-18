#include "compressiblecartesianpatch.h"
#include "ausmtools.h"

#define NI 100
#define NJ 10
#define NK 10

#include "code_blocks/minmod.h"
#include "code_blocks/cartesian_projection.h"

class TestPatch : public CompressibleCartesianPatch
{

public:

  virtual void computeResiduals();
  void correctBoundaries();

};

inline void TestPatch::computeResiduals()
{
  #include "code_blocks/compressible_residuals.h"
  #include "code_blocks/compressible_variables.h"

  for (size_t i = 0; i < NI - 1; ++i) {
    for (size_t j = 0; j < NJ; ++j) {
      for (size_t k = 0; k < NK; ++k) {

        #include "code_blocks/ausm_plus_x.h"

      }
    }
  }
}

inline void TestPatch::correctBoundaries()
{
  #include "code_blocks/compressible_variables.h"

  for (size_t j = 0; j < NJ; ++j) {
    for (size_t k = 0; k < NK; ++k) {
      setF(r, 0, j, k, f(r, 1, j, k));
      setF(r, sizeI() - 1, j, k, f(r, sizeI() - 2, j, k));
      setF(ru, 0, j, k, f(ru, 1, j, k));
      setF(ru, sizeI() - 1, j, k, f(ru, sizeI() - 2, j, k));
      setF(rv, 0, j, k, f(rv, 1, j, k));
      setF(rv, sizeI() - 1, j, k, f(rv, sizeI() - 2, j, k));
      setF(rw, 0, j, k, f(rw, 1, j, k));
      setF(rw, sizeI() - 1, j, k, f(rw, sizeI() - 2, j, k));
      setF(rE, 0, j, k, f(rE, 1, j, k));
      setF(rE, sizeI() - 1, j, k, f(rE, sizeI() - 2, j, k));
    }
  }
}


int main()
{
  TestPatch P;
  P.setupAligned(0, 0, 0, 1.0, 0.1, 0.1);
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
  while (t < 1e-4) {
    P.setOldState();
    for (int i_rk = 0; i_rk < 3; ++i_rk) {
      P.subStep(dt*alpha[i_rk]);
      P.correctBoundaries();
    }
    t += dt;
    cout << t << endl;
    break;
  }
  P.writeToVtk("test");
}

