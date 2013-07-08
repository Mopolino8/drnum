#ifndef SECONDORDER_H
#define SECONDORDER_H

#include "blockcfd.h"

struct SecondOrder
{
  static CUDA_DH real lim(real, real)
  {
    return 1;
  }
};

#endif // SECONDORDER_H
