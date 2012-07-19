#ifndef COMPRESSIBLECARTESIANPATCH_H
#define COMPRESSIBLECARTESIANPATCH_H

#include "cartesianpatch.h"

class CompressibleCartesianPatch : public CartesianPatch
{

public: // static attributes

  static const size_t i_new = 0;
  static const size_t i_old = 1;
  static const size_t i_res = 2;

public: // methods

  CompressibleCartesianPatch();

  void setState(size_t i, size_t j, size_t k, real p, real T);
  void setOldState();

  virtual void subStep(real dt);

};

#define COMPR_NEW_VARIABLES              \
  restrict real* r  = getVariable(i_new, 0); \
  restrict real* ru = getVariable(i_new, 1); \
  restrict real* rv = getVariable(i_new, 2); \
  restrict real* rw = getVariable(i_new, 3); \
  restrict real* rE = getVariable(i_new, 4); \

#define COMPR_OLD_VARIABLES                  \
  restrict real* old_r  = getVariable(i_old, 0); \
  restrict real* old_ru = getVariable(i_old, 1); \
  restrict real* old_rv = getVariable(i_old, 2); \
  restrict real* old_rw = getVariable(i_old, 3); \
  restrict real* old_rE = getVariable(i_old, 4); \

#define COMPR_RES_VARIABLES                  \
  restrict real* res_r  = getVariable(i_res, 0); \
  restrict real* res_ru = getVariable(i_res, 1); \
  restrict real* res_rv = getVariable(i_res, 2); \
  restrict real* res_rw = getVariable(i_res, 3); \
  restrict real* res_rE = getVariable(i_res, 4); \

#define CHECK_COMPR_FLUX {        \
  if (isnan(flux_r)) BUG;         \
  if (isinf(flux_r)) BUG;         \
  if (isnan(flux_ru)) BUG;        \
  if (isinf(flux_ru)) BUG;        \
  if (isnan(flux_rv)) BUG;        \
  if (isinf(flux_rv)) BUG;        \
  if (isnan(flux_rw)) BUG;        \
  if (isinf(flux_rw)) BUG;        \
  if (isnan(flux_rE)) BUG;        \
  if (isinf(flux_rE)) BUG;        \
}

#define ADD_COMPR_XFLUX             \
  f(res_r,  i, j, k)  -= flux_r;    \
  f(res_ru,  i, j, k) -= flux_ru;   \
  f(res_rv,  i, j, k) -= flux_rv;   \
  f(res_rw,  i, j, k) -= flux_rw;   \
  f(res_rE,  i, j, k) -= flux_rE;   \
  f(res_r,  i+1, j, k)  += flux_r;  \
  f(res_ru,  i+1, j, k) += flux_ru; \
  f(res_rv,  i+1, j, k) += flux_rv; \
  f(res_rw,  i+1, j, k) += flux_rw; \
  f(res_rE,  i+1, j, k) += flux_rE;

#endif // COMPRESSIBLECARTESIANPATCH_H
