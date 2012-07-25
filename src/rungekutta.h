#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "timeintegration.h"
#include <list>

class RungeKutta : public TimeIntegration
{

private: // attributes

  list<real> m_Alpha;


public:

  RungeKutta();

  void addAlpha(real a) { m_Alpha.push_back(a); }
  virtual void operator()(real dt);

};

#endif // RUNGEKUTTA_H
