#include "compressiblecartesianpatch.h"
#include "ausmtools.h"

#define NI 1000
#define NJ 2
#define NK 2

#include "reconstruction/upwind1.h"
#include "fluxes/ausmplus.h"

typedef CompressibleCartesianPatch<AusmPlus<Upwind1> > TestPatch;

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
  int sub_count = 0;
  write(P, count);
  while (t < 1e-3) {
    P.copyField(P.i_new, P.i_old);
    for (int i_rk = 0; i_rk < 3; ++i_rk) {
      P.subStep(dt*alpha[i_rk]);
    }
    t += dt;
    ++sub_count;
    if (sub_count == 10) {
      ++count;
      write(P, count);
      sub_count = 0;
    }
    real max_norm, l2_norm;
    P.computeVariableDifference(P.i_new, 0, P.i_old, 0, max_norm, l2_norm);
    cout << t << "  max: " << max_norm << "  L2: " << l2_norm << endl;
  }
}

