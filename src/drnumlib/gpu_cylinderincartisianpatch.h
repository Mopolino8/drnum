// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// +                                                                      +
// + This file is part of DrNUM.                                          +
// +                                                                      +
// + Copyright 2013 numrax GmbH, enGits GmbH                              +
// +                                                                      +
// + DrNUM is free software: you can redistribute it and/or modify        +
// + it under the terms of the GNU General Public License as published by +
// + the Free Software Foundation, either version 3 of the License, or    +
// + (at your option) any later version.                                  +
// +                                                                      +
// + DrNUM is distributed in the hope that it will be useful,             +
// + but WITHOUT ANY WARRANTY; without even the implied warranty of       +
// + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        +
// + GNU General Public License for more details.                         +
// +                                                                      +
// + You should have received a copy of the GNU General Public License    +
// + along with DrNUM. If not, see <http://www.gnu.org/licenses/>.        +
// +                                                                      +
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#ifndef CYLINDERINCARTISIANPATCH_H
#define CYLINDERINCARTISIANPATCH_H

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "perfectgas.h"

#include "drnum.h"
#include "genericoperation.h"
#include "gpu_cartesianpatch.h"

class GPU_CylinderInCartisianPatch;

class GPU_CylinderInCartisianPatch : public GenericOperation
{
  GPU_CartesianPatch *m_Patch;

  real    m_Tinf;
  real    m_Pinf;
  real    m_Rad;
  real    m_Rad2;
  real    m_Omega;
  vec3_t  m_Xo;
  vec3_t  m_Po;

  real   dot(vec3_t x1 , vec3_t x2);
  real   norm(vec3_t x);
  real   norm(real x, real y, real z);
  real   radiusSquare(vec3_t x1);
  real   scalarProduct(vec3_t x1 , vec3_t x2);
  vec3_t normalize(vec3_t x);
  vec3_t normalize(real x, real y, real z);

public:

  /**
   * This class sets two layers outside a cylinder of diameter m_Rad to a given pressure and velocity.
   * For storage during the computation of the layer data the residual field is abused.
   * The 5 storage indices per point are used as followed:
   *  - 1 cell type: 3-> in cylinder, 2-> layer 1, 1-> layer 2, 0-> outside, set to zero.
   *  - 2 total weight, w1 + w2 + w3.., accounting for the contribution of the cylinder edge cells to the "bc"
   *  - 3 weighted pressure, w1*p1 + w2*p2..
   *  - 4 u
   *  - 5 v
   */

  GPU_CylinderInCartisianPatch(GPU_CartesianPatch* patch);

  virtual void operator()();

};

GPU_CylinderInCartisianPatch::GPU_CylinderInCartisianPatch(GPU_CartesianPatch *patch)
{
  m_Patch = patch;
  m_Tinf  = 300.0;
  m_Pinf  = 1e5;
  m_Xo    = vec3_t(0.0, 0.0, 0.0);
  m_Rad   = 2.0;
  m_Rad2  = m_Rad*m_Rad;
  //m_Omega  = -1.04719775;
  m_Omega  = 17.3594;
  m_Po[0] = 0.5*m_Patch->dx();
  m_Po[1] = 0.5*m_Patch->dy();
  m_Po[2] = 0.5*m_Patch->dz();
  m_Po = m_Patch->getTransformInertial2This().transformReverse(m_Po);
}

inline real GPU_CylinderInCartisianPatch::norm(vec3_t x)
{
  return sqrt(x[0]*x[0] + x[1]*x[1] +x[2]*x[2]);
}

inline real GPU_CylinderInCartisianPatch::norm(real x, real y, real z)
{
  return sqrt(x*x + y*y + z*z);
}

inline vec3_t GPU_CylinderInCartisianPatch::normalize(real x, real y, real z)
{
  vec3_t x_return;
  real norm = norm(x, y, z);
  x_return[0] = x/norm;
  x_return[1] = y/norm;
  x_return[2] = z/norm;

  return x_return;
}

inline vec3_t GPU_CylinderInCartisianPatch::normalize(vec3_t x)
{
  vec3_t x_return;
  x_return[0] = x[0]/norm(x);
  x_return[1] = x[1]/norm(x);
  x_return[2] = x[2]/norm(x);

  return  x;
}

inline real GPU_CylinderInCartisianPatch::radiusSquare(vec3_t x1)
{
  return GPU_CylinderInCartisianPatch::radiusSquare(x1[0], x1[1], x1[2]);
}

inline real GPU_CylinderInCartisianPatch::radiusSquare(real x, real y, real z)
{
  return (x*x + y*y + z*z);
}

inline real GPU_CylinderInCartisianPatch::dot(vec3_t x1, vec3_t x2)
{
  return (x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2]);
}

inline real GPU_CylinderInCartisianPatch::scalarProduct(vec3_t x1, vec3_t x2)
{
  return GPU_CylinderInCartisianPatch::dot(x1, x2)/GPU_CylinderInCartisianPatch::norm(x2);
}

template <unsigned int DIM>
inline void GPU_CylinderInCartisianPatch::computeDpAndWeight(GPU_CartesianPatch* patch, vec3_t x_PmO, vec3_t x_P, size_t i, size_t j, size_t k, real dx, real& dp, real& weight) {
  vec3_t x_im;
  x_im[0] = (real(i) + 0.5)*patch.dx() + x_PmO[0];
  x_im[1] = (real(j) + 0.5)*patch.dy() + x_PmO[1];
  x_im[2] = (real(k) + 0.5)*patch.dz() + x_PmO[2];

  real var1[5];
  real T, p, u, v, w;
  patch.getVar(DIM, 0, i, j, k, var1);
  PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);

  weight = scalarProduct(normalize(x_im), normalize(x_p));
  // dp = dr*v^2*rho/r
  // Here dr = dx*weight
  dp     = weight*p + weight*weight*m_Patch.dx()*(u*u + v*v + w*w)*var1[0]/norm(x_im);
}


template <unsigned int DIM>
__global__ void GPU_CylinderInCartisianPatch_Mark(GPU_CartesianPatch patch, vec3_t x_PmO, unsigned int radius2)
{
  size_t i = blockDim.y*blockIdx.x + threadIdx.y;
  size_t j = 2*blockIdx.y;
  size_t k = threadIdx.x;
  if (j >= patch.sizeJ() || i >= patch.sizeI()) return;

  //real x = (real)0.5*patch.dx() + i*patch.dx();
  //real y = (real)0.5*patch.dy() + j*patch.dy();
  //real z = (real)0.5*patch.dz() + k*patch.dz();
  real x = (real(i) + 0.5)*patch.dx() + x_PmO[0];
  real y = (real(j) + 0.5)*patch.dy() + x_PmO[1];
  real z = (real(k) + 0.5)*patch.dz() + x_PmO[2];

  if (radiusSquare(x, y, z) < radius2) {
    patch.f(2, 0, i, j, k) = 3;
  }
}

template <unsigned int DIM>
__global__ void GPU_CylinderInCartisianPatch_ComputeLayer(GPU_CartesianPatch patch, vec3_t x_PmO, unsigned int radius2, unsigned int layer, real omega)
{
  size_t i = 1 + blockDim.y*blockIdx.x + threadIdx.y;
  size_t j = 1 + 2*blockIdx.y;
  size_t k = threadIdx.x;
  if (j >= (patch.sizeJ() - 1) || i >= (patch.sizeI() - 1)) return;

  real l_p = patch.f(2, 0, i, j, k);
  // if current point == 3, ie in cylinder, return
  if ( l_p > 2.5 ) return;
  vec3_t x_p;
  x_p[0] = (real(i) + 0.5)*patch.dx() + x_PmO[0];
  x_p[1] = (real(j) + 0.5)*patch.dy() + x_PmO[1];
  x_p[2] = (real(k) + 0.5)*patch.dz() + x_PmO[2];

  real l_im = patch.f(2, 0, i - 1, j    , k);
  real l_ip = patch.f(2, 0, i + 1, j    , k);
  real l_jm = patch.f(2, 0, i    , j - 1, k);
  real l_jp = patch.f(2, 0, i    , j + 1, k);

  real T, p, u, v, w;
  real var1[5];
  real p_sum = 0, w_sum = 0;
  bool updated = false;

  // if l_im point in cylinder, interpolate to l_p
  if (l_im > (3.5 - layer) && l_p < 0.5) {
    vec3_t x_im;
    x_im[0] = (real(i - 1) + 0.5)*patch.dx() + x_PmO[0];
    x_im[1] = (real(j)     + 0.5)*patch.dy() + x_PmO[1];
    x_im[2] = (real(k)     + 0.5)*patch.dz() + x_PmO[2];

    patch.getVar(DIM, 2, i - 1, j, k, var1);
    PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);

    real weight = scalarProduct(normalize(x_im), normalize(x_p));
    // dp = dr*v^2*rho/r
    // Here dr = dx*weight
    real dp     = weight*p + weight*weight*m_Patch.dx()*(u*u + v*v + w*w)*var1[0]/norm(x_im);
    p_sum += dp;
    w_sum += weight;

    updated = true;
  }
  // if l_ip point in cylinder, interpolate to l_p
  if (l_ip > (3.5 - layer) && l_p < 0.5) {
    vec3_t x_ip;
    x_ip[0] = (real(i + 1) + 0.5)*patch.dx() + x_PmO[0];
    x_ip[1] = (real(j)     + 0.5)*patch.dy() + x_PmO[1];
    x_ip[2] = (real(k)     + 0.5)*patch.dz() + x_PmO[2];

    patch.getVar(DIM, 2, i + 1, j, k, var1);
    PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);

    real weight = scalarProduct(normalize(x_ip), normalize(x_p));
    // dp = dr*v^2*rho/r
    // Here dr = dx*weight
    real dp     = weight*p + weight*weight*m_Patch.dx()*(u*u + v*v + w*w)*var1[0]/norm(x_ip);
    p_sum += dp;
    w_sum += weight;

    updated = true;
  }
  // if l_jm point in cylinder, interpolate to l_p
  if (l_jm > (3.5 - layer) && l_p < 0.5) {
    vec3_t x_jm;
    x_jm[0] = (real(i)     + 0.5)*patch.dx() + x_PmO[0];
    x_jm[1] = (real(j - 1) + 0.5)*patch.dy() + x_PmO[1];
    x_jm[2] = (real(k)     + 0.5)*patch.dz() + x_PmO[2];

    patch.getVar(DIM, 2, i, j - 1, k, var1);
    PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);

    real weight = scalarProduct(normalize(x_jm), normalize(x_p));
    // dp = dr*v^2*rho/r
    // Here dr = dy*weight
    real dp     = weight*p + weight*weight*m_Patch.dy()*(u*u + v*v + w*w)*var1[0]/norm(x_jm);
    p_sum += dp;
    w_sum += weight;

    updated = true;
  }
  // if l_jp point in cylinder, interpolate to l_p
  if (l_ip > (3.5 - layer) && l_p < 0.5) {
    vec3_t x_jp;
    x_jp[0] = (real(i)     + 0.5)*patch.dx() + x_PmO[0];
    x_jp[1] = (real(j + 1) + 0.5)*patch.dy() + x_PmO[1];
    x_jp[2] = (real(k)     + 0.5)*patch.dz() + x_PmO[2];

    patch.getVar(DIM, 2, i, j + 1, k, var1);
    PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);

    real weight = scalarProduct(normalize(x_jp), normalize(x_p));
    // dp = dr*v^2*rho/r
    // Here dr = dy*weight
    real dp     = weight*p + weight*weight*m_Patch.dy()*(u*u + v*v + w*w)*var1[0]/norm(x_jp);
    p_sum += dp;
    w_sum += weight;

    updated = true;
  }

  if (updated == true) {
    patch.f(2, 0, i, j, k) = 3 - layer;

    real u_layer =  omega*x_p[1];         // ==  m_Omega*norm(r_2)*r_2[1]/norm(r_2) ->  v_t*r*sin(th);
    real v_layer = -omega*x_p[0];         // == -m_Omega*norm(r_2)*r_2[0]/norm(r_2) -> -v_t*r*cos(th);

    PerfectGas::primitiveToConservative(p_sum/w_sum, m_Tinf, u_layer, v_layer, 0, d_var1);
    m_Patch->setVar(dim, 0, i, j, k, d_var1);
  }
  else if (patch.f(2, 0, i, j, k) < 0.5) {
    PerfectGas::primitiveToConservative(m_Pinf, m_Tinf, 0, 0, 0, d_var1);
    m_Patch->setVar(dim, 0, i, j, k, d_var1);
  }
}

void GPU_CylinderInCartisianPatch::operator ()()
{
  dim_t<5> dim;

  vec3_t x_PmO = m_Po - m_Xo;

  real var1[5];
  for (int i_var = 0; i_var < 5; ++i_var) {
    var1[i_var] = 0;
  }
  // Clear residual field
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        m_Patch->setVar(dim, 2, i, j, k, var1);
      }
    }
  }



  vec3_t r, x_p;
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        x_p[0] = (real(i) + 0.5)*m_Patch->dx();
        x_p[1] = (real(j) + 0.5)*m_Patch->dy();
        x_p[2] = (real(k) + 0.5)*m_Patch->dz();
        x_p = m_Patch->getTransformInertial2This().transformReverse(x_p);

        r[0] = x_p[0] - m_XOrg[0];
        r[1] = x_p[1] - m_XOrg[1];
        r[2] = x_p[2] - m_XOrg[2];

        if ((r[0]*r[0] + r[1]*r[1]) < m_Rad*m_Rad) {
          // set to interior
          var1[0] = 1;
          // Abuse the residual field for storage
          m_Patch->setVar(dim, 2, i, j, k, var1);
          var1[0] = 0;
        }
        else {
          // clear data field for cells outside of cylinder
          m_Patch->setVar(dim, 0, i, j, k, var1);
        }
      }
    }
  }

  //
  real var2[5], d_var1[5], d_var2[5];
  vec3_t r_1, x_p1, r_2, x_p2;
  for (int i_var = 0; i_var < 5; ++i_var) {
    var1[i_var]   = 0;  // Residual field, pt 1
    var2[i_var]   = 0;  // Residual field, pt 2
    d_var1[i_var] = 0;  // Data field, pt 1
    d_var2[i_var] = 0;  // Data field, pt 2
  }
  for (size_t layer = 0; layer < 2; ++layer) {
    // iterate over I values, and compare i, i-1
    for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
      for (size_t i = 1; i < m_Patch->sizeI(); ++i) {
        for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
          m_Patch->getVar(dim, 2, i - 1, j, k, var1);
          m_Patch->getVar(dim, 2, i    , j, k, var2);

          // 0.5 to ensure real / int comparison works
          if ((var1[0] > (0.5 + layer) && var2[0] <  0.5 ) ||
              (var1[0] <  0.5          && var2[0] > (0.5 + layer)))
          {
            x_p1[0] = (real(i - 1) + 0.5)*m_Patch->dx();
            x_p1[1] = (real(j)     + 0.5)*m_Patch->dy();
            x_p1[2] = (real(k)     + 0.5)*m_Patch->dz();
            x_p1 = m_Patch->getTransformInertial2This().transformReverse(x_p1);
            r_1[0] = x_p1[0] - m_XOrg[0];
            r_1[1] = x_p1[1] - m_XOrg[1];
            r_1[2] = x_p1[2] - m_XOrg[2];

            x_p2[0] = (real(i) + 0.5)*m_Patch->dx();
            x_p2[1] = (real(j) + 0.5)*m_Patch->dy();
            x_p2[2] = (real(k) + 0.5)*m_Patch->dz();
            x_p2 = m_Patch->getTransformInertial2This().transformReverse(x_p2);
            r_2[0] = x_p2[0] - m_XOrg[0];
            r_2[1] = x_p2[1] - m_XOrg[1];
            r_2[2] = x_p2[2] - m_XOrg[2];

            m_Patch->getVar(dim, 0, i - 1, j, k, d_var1);
            m_Patch->getVar(dim, 0, i    , j, k, d_var2);

            if (var1[0] > (0.5 + layer)) {
              // pnt 1 is reference, interpolate to pnt 2
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var1, p, T, u, v, w);

              real u_layer =  m_Omega*r_2[1];         // ==  m_Omega*norm(r_2)*r_2[1]/norm(r_2) ->  v_t*r*sin(th);
              real v_layer = -m_Omega*r_2[0];         // == -m_Omega*norm(r_2)*r_2[0]/norm(r_2) -> -v_t*r*cos(th);
              real weight  = scalarProduct(normalize(r_1), normalize(r_2));
              real dp_new  = weight*(p + m_Patch->dx()*(u*u + v*v + w*w)*d_var1[0]/norm(r_1));

              var2[0]  = -1;     // mark layer
              var2[1] += weight;
              var2[2] += dp_new;
              var2[3]  = u_layer;
              var2[4]  = v_layer;
              m_Patch->setVar(dim, 2, i, j, k, var2);
            } else {
              // pnt 2 is reference, extrapolate to pnt 1
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var2, p, T, u, v, w);

              real u_layer =  m_Omega*r_1[1];         // == see above
              real v_layer = -m_Omega*r_1[0];         // == see above
              real weight  = scalarProduct(normalize(r_2), normalize(r_1));
              real dp_new  = weight*(p + m_Patch->dx()*(u*u + v*v + w*w)*d_var2[0]/norm(r_2));

              var1[0]  = -1;     // mark layer
              var1[1] += weight;
              var1[2] += dp_new;
              var1[3]  = u_layer;
              var1[4]  = v_layer;
              m_Patch->setVar(dim, 2, i - 1, j, k, var1);
            }
          } // end if loop
        } // end i loop
      } // end j loop
    } // end k loop

    // J values
    for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
      for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
        for (size_t j = 1; j < m_Patch->sizeJ(); ++j) {
          m_Patch->getVar(dim, 2, i, j - 1, k, var1);
          m_Patch->getVar(dim, 2, i, j    , k, var2);

          if ((var1[0] > (0.5 + layer) && var2[0] <  0.5) ||
              (var1[0] <  0.5          && var2[0] > (0.5 + layer)))
          {
            x_p1[0] = (real(i)     + 0.5)*m_Patch->dx();
            x_p1[1] = (real(j - 1) + 0.5)*m_Patch->dy();
            x_p1[2] = (real(k)     + 0.5)*m_Patch->dz();
            x_p1 = m_Patch->getTransformInertial2This().transformReverse(x_p1);
            r_1[0] = x_p1[0] - m_XOrg[0];
            r_1[1] = x_p1[1] - m_XOrg[1];
            r_1[2] = x_p1[2] - m_XOrg[2];

            x_p2[0] = (real(i) + 0.5)*m_Patch->dx();
            x_p2[1] = (real(j) + 0.5)*m_Patch->dy();
            x_p2[2] = (real(k) + 0.5)*m_Patch->dz();
            x_p2 = m_Patch->getTransformInertial2This().transformReverse(x_p2);
            r_2[0] = x_p2[0] - m_XOrg[0];
            r_2[1] = x_p2[1] - m_XOrg[1];
            r_2[2] = x_p2[2] - m_XOrg[2];

            m_Patch->getVar(dim, 0, i, j - 1, k, d_var1);
            m_Patch->getVar(dim, 0, i, j    , k, d_var2);

            if (var1[0] > (0.5 + layer) && var2[0] < (0.5 + layer)) {
              // pnt 1 is reference, extrapolate to pnt 2
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var1, p, T, u, v, w);

              real u_layer =  m_Omega*r_2[1];         // ==  m_Omega*norm(r_2)*r_2[1]/norm(r_2) ->  v_t*r*sin(th);
              real v_layer = -m_Omega*r_2[0];         // == -m_Omega*norm(r_2)*r_2[0]/norm(r_2) -> -v_t*r*cos(th);
              real weight  = scalarProduct(normalize(r_1), normalize(r_2));
              real dp_new  = weight*(p + m_Patch->dy()*(u*u + v*v + w*w)*d_var1[0]/norm(r_1));

              var2[0]  = -1;     // mark layer
              var2[1] += weight;
              var2[2] += dp_new;
              var2[3]  = u_layer;
              var2[4]  = v_layer;
              m_Patch->setVar(dim, 2, i, j, k, var2);
            }
            else if (var1[0] < (0.5 + layer) && var2[0] > (0.5 + layer)) {
              // pnt 2 is reference, extrapolate to pnt 1
              real T, p, u, v, w;
              PerfectGas::conservativeToPrimitive(d_var2, p, T, u, v, w);

              real u_layer =  m_Omega*r_1[1];         // == see above
              real v_layer = -m_Omega*r_1[0];         // == see above
              real weight = scalarProduct(normalize(r_2), normalize(r_1));
              real dp_new = weight*(p + m_Patch->dy()*(u*u + v*v + w*w)*d_var2[0]/norm(r_2));

              var1[0] = -1;     // mark layer
              var1[1] += weight;
              var1[2] += dp_new;
              var1[3] = u_layer;
              var1[4] = v_layer;
              m_Patch->setVar(dim, 2, i, j - 1, k, var1);
            }
          } // end if loop
        } // end i loop
      } // end j loop
    } // end k loop

    // Set the marked layers to the correct type number. 2 for layer 1, 3 for layer 2.
    // Then update the data field
    for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
      for (size_t j = 1; j < m_Patch->sizeJ(); ++j) {
        for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
          m_Patch->getVar(dim, 2, i, j, k, var1);
          if (var1[0] < -0.5) {
            var1[0] = 2 + layer;
            m_Patch->setVar(dim, 2, i, j, k, var1);

            m_Patch->getVar(dim, 0, i, j, k, d_var1);
            real u, v, p, weight;
            weight = var1[1];
            p      = var1[2];
            u      = var1[3];
            v      = var1[4];
            PerfectGas::primitiveToConservative(p/weight, m_Temp, u, v, 0, d_var1);
            m_Patch->setVar(dim, 0, i, j, k, d_var1);
          }
        } // end i loop
      } // end j loop
    } // end k loop

  } // end layer loop


  // Set the outside field to zero.
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        m_Patch->getVar(dim, 2, i, j, k, var1);
        x_p2[0] = (real(i) + 0.5)*m_Patch->dx();
        x_p2[1] = (real(j) + 0.5)*m_Patch->dy();
        x_p2[2] = (real(k) + 0.5)*m_Patch->dz();
        x_p2 = m_Patch->getTransformInertial2This().transformReverse(x_p2);
        if (var1[0] < 0.5) {
          PerfectGas::primitiveToConservative(1e5, m_Temp, 0, 0, 0, d_var1);
          m_Patch->setVar(dim, 0, i, j, k, d_var1);
        }
        if (var1[0] < 1.5 && var1[0] > 0.5) {
          cout << "data1 " << x_p2[0]
               << " " << x_p2[1]
               << " " << x_p2[2]
               << endl;
        }
        if (var1[0] < 2.5 && var1[0] > 1.5) {
          cout << "data2 " << x_p2[0]
               << " " << x_p2[1]
               << " " << x_p2[2]
               << endl;
        }
        if (var1[0] > 2.5) {
          cout << "data3 " << x_p2[0]
               << " " << x_p2[1]
               << " " << x_p2[2]
               << endl;
        }
      } // end i loop
    } // end j loop
  } // end k loop

  // Clear residual field
  for (int i_var = 0; i_var < 5; ++i_var) {
    var1[i_var] = 0;
  }
  for (size_t i = 0; i < m_Patch->sizeI(); ++i) {
    for (size_t j = 0; j < m_Patch->sizeJ(); ++j) {
      for (size_t k = 0; k < m_Patch->sizeK(); ++k) {
        m_Patch->setVar(dim, 2, i, j, k, var1);
      }
    }
  }

  m_Patch->copyFieldToDevice(0);
}

#endif // CYLINDERINCARTISIANPATCH_H
