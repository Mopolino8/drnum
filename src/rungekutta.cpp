#include "rungekutta.h"

RungeKutta::RungeKutta()
{
}

void RungeKutta::operator()(real dt)
{
  bool first_step = true;
  copyDonorData(0);
  copyField(0, 1);
  for (list<real>::iterator i = m_Alpha.begin(); i != m_Alpha.end(); ++i) {
    if (!first_step) {
      copyDonorData(0);
    }
    first_step = false;
    computeIterators((*i)*dt);
    runPostOperations();
    countFlops(1);
  }
}
