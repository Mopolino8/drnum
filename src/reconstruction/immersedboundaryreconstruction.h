#ifndef IMMERSEDBOUNDARYRECONSTRUCTION_H
#define IMMERSEDBOUNDARYRECONSTRUCTION_H

#include "blockcfd.h"
#include  "perfectgas.h"

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

    // auxilary variable sets
    real* var0 = new real [num_vars];
    real* var1 = new real [num_vars];
    real* bvar = new real [num_vars];

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

    countFlops(9);

    // compute corrected values on the surface of the shape
    weight *= 2;
    m_Reconstruction.project(patch, var, i_field, num_vars, i1, j1, k1, i2, j2, k2);

    // check if it is a boundary flux
    // stop here if yes and if the intersection is outside of the domain (patch)
    // ... ==> default to normal reconstruction
    //
    if (!patch->checkRange(i2, j2, k2) && weight > 1) {
      return;
    }

    patch->getVar(i_field, i1, j1, k1, var1);
    for (size_t i_var = 0; i_var < num_vars; ++i_var) {
      bvar[i_var] = weight*var[i_var] + (1-weight)*var1[i_var];
    }
    TBoundaryCondition::correct(nx, ny, nz, bvar);
    countFlops(1 + 4*num_vars);

    // compute the indices of one node further away from the boundary surface
    size_t i0 = 2*i1 - i2;
    size_t j0 = 2*j1 - j2;
    size_t k0 = 2*k1 - k2;
    if (!patch->checkRange(i0, j0, k0)) {
      BUG;
    }

    // extrapolate (interpolate) to the face i + 1/2
    patch->getVar(i_field, i0, j0, k0, var0);
    for (size_t i_var = 0; i_var < num_vars; ++i_var) {
      //var[i_var] = final_weight*bvar[i_var] + (1 - final_weight)*var0[i_var];
      var[i_var] = bvar[i_var] + 0.5*(1 - weight)*(bvar[i_var] - var0[i_var])/(1 + 0.5*weight);
      //var[i_var] = bvar[i_var];// + (1 - weight)*(bvar[i_var] - 0.5*(var0[i_var] + var1[i_var])/(1 + weight));
    }

    countFlops(8 + 4*num_vars);

    delete [] var0;
    delete [] var1;
    delete [] bvar;
    return;
  }

  if (m_Shape->isInside(2*x1 - x2, 2*y1 - y2, 2*z1 - z2)) {
    /// @todo look at improving this
    /*
    if (!patch->checkRange(i2, j2, k2)) {
      // first order upwind ...
      patch->getVar(i_field, i1, j1, k1, var);
    } else {
      // central ...
      real* var1 = new real [num_vars];
      real* var2 = new real [num_vars];
      patch->getVar(i_field, i1, j1, k1, var1);
      patch->getVar(i_field, i2, j2, k2, var2);
      for (size_t i_var = 0; i_var < num_vars; ++i_var) {
        var[i_var] = 0.5*(var1[i_var] + var2[i_var]);
      }
    }
    */
    // first order upwind ...
    patch->getVar(i_field, i1, j1, k1, var);

  } else {
    // default to normal reconstruction
    m_Reconstruction.project(patch, var, i_field, num_vars, i1, j1, k1, i2, j2, k2);
  }
}

#endif // IMMERSEDBOUNDARYRECONSTRUCTION_H
