#include "rungekutta.h"

RungeKutta::RungeKutta()
{
}

void RungeKutta::operator()(real dt)
{
  copyField(0, 1);
  for (list<real>::iterator i = m_Alpha.begin(); i != m_Alpha.end(); ++i) {
    computeIterators((*i)*dt);
    countFlops(1);
  }
}
