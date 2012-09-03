#ifndef SECONDORDER_H
#define SECONDORDER_H

#include "blockcfd.h"

struct SecondOrder
{
  static real delta(real delta1, real)
  {
    return delta1;
  }
};

#endif // SECONDORDER_H
