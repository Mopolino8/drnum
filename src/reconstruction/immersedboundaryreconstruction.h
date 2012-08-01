#ifndef IMMERSEDBOUNDARYRECONSTRUCTION_H
#define IMMERSEDBOUNDARYRECONSTRUCTION_H

#include "blockcfd.h"

template <typename TReconstruction, typename TShape, typename TBoundaryCondition>
class ImmersedBoundaryReconstruction
{

  TShape* m_Shape;
  TReconstruction m_Reconstruction;


public:

  ImmersedBoundaryReconstruction(TShape* shape) { m_Shape = shape; }

  void project(CartesianPatch *patch, real* var, size_t i_field, size_t num_vars,
               size_t i1, size_t j1, size_t k1,
               size_t i2, size_t j2, size_t k2,
               real x1, real y1, real z1,
               real x2, real y2, real z2);
};

template <typename TReconstruction, typename TShape, typename TBoundaryCondition>
inline void ImmersedBoundaryReconstruction<TReconstruction, TShape, TBoundaryCondition>::project
(
  CartesianPatch *patch, real *var, size_t i_field, size_t num_vars,
  size_t i1, size_t j1, size_t k1,
  size_t i2, size_t j2, size_t k2,
  real x1, real y1, real z1,
  real x2, real y2, real z2
)
{
  if (patch->checkRange(i1, j1, k1) && patch->checkRange(i2, j2, k2)) {
    real k, nx, ny, nz;

    if (m_Shape->getBoundaryMetric(x1, y1, z1, x2, y2, z2, k, nx, ny, nz))
    {
      // auxilary variable sets
      real* var0 = new real [num_vars];
      real* var1 = new real [num_vars];
      real* bvar = new real [num_vars];

      // make sure we always extrapolate from outside the shape
      if ((x2-x1)*nx + (y2-y1)*ny + (z2-z1)*nz > 0) {
        swap(i1, i2);
        swap(j1, j2);
        swap(k1, k2);
        k = 1 - k;
      }
      countFlops(9);

      // compute corrected values on the surface of the shape
      k *= 2;
      m_Reconstruction.project(patch, var, i_field, num_vars, i1, j1, k1, i2, j2, k2);
      patch->getVar(i_field, i1, j1, k1, var1);
      for (size_t i_var = 0; i_var < num_vars; ++i_var) {
        bvar[i_var] = k*var[i_var] + (1-k)*var1[i_var];
      }
      TBoundaryCondition::correct(nx, ny, nz, bvar);
      countFlops(1 + 4*num_vars);

      // extrapolate (interpolate) to the face 1 + 1/2
      k = 3.0/(2 + k);
      size_t i0 = 2*i1 - i2;
      size_t j0 = 2*j1 - j2;
      size_t k0 = 2*k1 - k2;

      if (!patch->checkRange(i0, j0, k0)) {
        BUG;
      }

      patch->getVar(i_field, i0, j0, k0, var0);
      for (size_t i_var = 0; i_var < num_vars; ++i_var) {
        var[i_var] = k*bvar[i_var] + (1-k)*var0[i_var];
      }
      countFlops(8 + 4*num_vars);

      delete [] var0;
      delete [] var1;
      delete [] bvar;
      return;
    }
  }

  // default to normal reconstruction
  m_Reconstruction.project(patch, var, i_field, num_vars, i1, j1, k1, i2, j2, k2);
}

#endif // IMMERSEDBOUNDARYRECONSTRUCTION_H
