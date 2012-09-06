#ifndef SECONDORDER_H
#define SECONDORDER_H

#include "blockcfd.h"

struct SecondOrder
{
  static real lim(real delta1, real)
  {
    return 1;
  }
};

#endif // SECONDORDER_H
