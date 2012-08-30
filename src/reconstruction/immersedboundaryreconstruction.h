#ifndef IMMERSEDBOUNDARYRECONSTRUCTION_H
#define IMMERSEDBOUNDARYRECONSTRUCTION_H

#include "blockcfd.h"

/**
 * A genric reconstruction which respects immersed boundaries.
 * The boundary description has to be passed as a template parameter (TShape)
 */
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
  real weight, nx, ny, nz;

  // check if the edge (i1, j1, k1) -> (i2, j2, k2) crosses an immersed boundary
  if (m_Shape->getBoundaryMetric(x1, y1, z1, x2, y2, z2, weight, nx, ny, nz)) {

    // make sure we always extrapolate from outside the shape (i.e. from the fluid part)
    if ((x2-x1)*nx + (y2-y1)*ny + (z2-z1)*nz > 0) {
      swap(i1, i2);
      swap(j1, j2);
      swap(k1, k2);
      weight = 1 - weight;
    }
    if (!patch->checkRange(i1, j1, k1)) {
      BUG;
    }
    patch->getVar(i_field, i1, j1, k1, var);
    TBoundaryCondition::correct(nx, ny, nz, var);
    return;
  }

  if (m_Shape->isInside(2*x1 - x2, 2*y1 - y2, 2*z1 - z2)) {

    // first order upwind ...
    patch->getVar(i_field, i1, j1, k1, var);

  } else {

    // default to normal reconstruction
    m_Reconstruction.project(patch, var, i_field, num_vars, i1, j1, k1, i2, j2, k2);

  }
}

#endif // IMMERSEDBOUNDARYRECONSTRUCTION_H
