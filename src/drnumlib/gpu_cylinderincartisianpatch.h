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
#ifndef GPU_CYLINDERINCARTISIANPATCH_H
#define GPU_CYLINDERINCARTISIANPATCH_H

#ifdef CUDA

#include <cuda.h>
#include <cuda_runtime_api.h>

#include "perfectgas.h"

#include "drnum.h"
#include "genericoperation.h"
#include "gpu_cartesianpatch.h"

class GPU_CylinderInCartisianPatch;

class GPU_CylinderInCartisianPatch : public GenericOperation
{
  GPU_CartesianPatch m_Patch;

  real   m_Tinf;
  real   m_Pinf;
  real   m_Rad;
  real   m_Omega;
  real   m_Height;
  vec3_t m_Xo;
  vec3_t m_Po;
  size_t m_MaxNumThreads;

public:

  /**
   * This class sets two layers outside a cylinder of diameter m_Rad to a given pressure and velocity.
   * For storage during the computation of the layer data the residual field is abused.
   * The 5 storage indices per point are used as followed:
   *  - index 1 : cell type => 3-> in cylinder, 2-> layer 1, 1-> layer 2, 0-> outside, set to zero.
   */

  GPU_CylinderInCartisianPatch(CartesianPatch *patch, int cuda_device, size_t thread_limit);

  virtual void operator()();

};

GPU_CylinderInCartisianPatch::GPU_CylinderInCartisianPatch(CartesianPatch *patch, int cuda_device, size_t thread_limit) : m_Patch(patch)
{
  m_Tinf   = 300.0;
  m_Pinf   = 1e5;
  m_Xo     = vec3_t(0.0, 0.0, 0.0);
  m_Rad    = 1.0;
  m_Height = 0.5;
  //m_Omega  = -1.04719775;
  m_Omega  = 17.3594;
  m_Po[0] = 0.5*patch->dx();
  m_Po[1] = 0.5*patch->dy();
  m_Po[2] = 0.5*patch->dz();
  m_Po    = patch->getTransformInertial2This().transformReverse(m_Po);

  int count;
  if (cudaGetDeviceCount(&count) != cudaSuccess) {
    cerr << "error detecting CUDA devices" << endl;
    exit(EXIT_FAILURE);
  }
  if (count < cuda_device + 1) {
    CudaTools::info();
    cerr << "specified CUDA device does not exists" << endl;
    exit(EXIT_FAILURE);
  }
  cudaDeviceProp prop;
  if (cudaGetDeviceProperties(&prop, cuda_device) != cudaSuccess) {
    cerr << "error fetching device properties" << endl;
    exit(EXIT_FAILURE);
  }
  cudaSetDevice(cuda_device);
  m_MaxNumThreads = min(prop.maxThreadsPerBlock, prop.maxThreadsDim[0]);
  if (thread_limit > 0) {
    m_MaxNumThreads = min(thread_limit, m_MaxNumThreads);
  }
  //cudaDeviceSetLimit( cudaLimitPrintfFifoSize, 32 * 1024 * 1024 );
}

CUDA_DH inline real vectorSquare(real x, real y, real z)
{
  return (x*x + y*y + z*z);
}

CUDA_DH inline real norm(real x, real y, real z)
{
  return sqrt(x*x + y*y + z*z);
}

CUDA_DH inline real dot(real x1, real y1, real z1, real x2, real y2, real z2)
{
  return (x1*x2 + x1*x2 + x1*x2);
}

CUDA_DH inline real scalarProduct(real x1, real y1, real z1, real x2, real y2, real z2)
{
  return dot(x1, y1, z1, x2, y2, z2)/norm(x2, y2, z2);
}

CUDA_DH inline real normalizedScalarProduct(real x1, real y1, real z1, real x2, real y2, real z2)
{
  real norm1 = norm(x1, y1, z1);
  x1 /= norm1;
  y1 /= norm1;
  z1 /= norm1;

  real norm2 = norm(x2, y2, z2);
  x2 /= norm2;
  y2 /= norm2;
  z2 /= norm2;

  return scalarProduct(x1, y1, z1, x2, y2, z2);
}

CUDA_DH inline void computeDpAndWeight(GPU_CartesianPatch& patch, real x_PmO, real y_PmO, real z_PmO, real x_p, real y_p, real z_p, size_t i, size_t j, size_t k, real dx, real& dp, real& weight)
{
  dim_t<5> dim;
  real x = (real(i) + 0.5)*patch.dx() + x_PmO;
  real y = (real(j) + 0.5)*patch.dy() + y_PmO;
  real z = (real(k) + 0.5)*patch.dz() + z_PmO;

  real var1[5];
  real T, p, u, v, w;
  patch.getVar(dim, 0, i, j, k, var1);
  PerfectGas::conservativeToPrimitive(var1, p, T, u, v, w);

  weight = normalizedScalarProduct(x, y, z, x_p, y_p, z_p);

  // dp = dr.v_th^2.rho/r
  // Here dr = dx*weight
  real v_n   = scalarProduct(u, v, 0, x_p, y_p, 0);
  real v_th2 = vectorSquare(u, v, 0) - v_n*v_n;
  //dp = weight*p + weight*weight*dx*(u*u + v*v + w*w)*var1[0]/norm(x_im);
  dp = weight*p + weight*weight*dx*v_th2*var1[0]/norm(x, y, z);
}

__global__ void GPU_CylinderInCartisianPatch_Mark(GPU_CartesianPatch patch, real x_PmO, real y_PmO, real z_PmO, real radius2, real height2)
{
  //size_t i = blockDim.y*blockIdx.x + threadIdx.y;
  //size_t j = 2*blockIdx.y;
  //size_t k = threadIdx.x;
  size_t i = blockIdx.x;
  size_t j = blockIdx.y;
  size_t k = threadIdx.x;

  if (i >= patch.sizeI() || j >= patch.sizeJ() || k >= patch.sizeK() ) return;

  real x = (real(i) + 0.5)*patch.dx() + x_PmO;
  real y = (real(j) + 0.5)*patch.dy() + y_PmO;
  real z = (real(k) + 0.5)*patch.dz() + z_PmO;

  //if (vectorSquare(x, y, 0.) < radius2) {
  if (vectorSquare(x, y, 0.) < radius2 && z*z < height2) {
    patch.f(2, 0, i, j, k) = 3;
  }
  else {
    patch.f(2, 0, i, j, k) = 0;
  }
}

__global__ void GPU_CylinderInCartisianPatch_ComputeLayer(GPU_CartesianPatch patch, real x_PmO, real y_PmO, real z_PmO, real p_inf, real t_inf, real omega, unsigned int layer)
{
  dim_t<5> dim;

  size_t i = 1 + blockIdx.x;
  size_t j = 1 + blockIdx.y;
  size_t k = 1 + threadIdx.x;

  if (i >= (patch.sizeI() - 1) || j >= (patch.sizeJ() - 1) || k >= (patch.sizeK() - 1)) return;

  real l_p = patch.f(2, 0, i, j, k);
  // if current point marked, ie in cylinder, return
  if ( l_p > 0.5 ) return;

  real l_im = patch.f(2, 0, i - 1, j    , k    );
  real l_ip = patch.f(2, 0, i + 1, j    , k    );
  real l_jm = patch.f(2, 0, i    , j - 1, k    );
  real l_jp = patch.f(2, 0, i    , j + 1, k    );
  real l_km = patch.f(2, 0, i    , j    , k - 1);
  real l_kp = patch.f(2, 0, i    , j    , k + 1);

  real x_p = (real(i) + 0.5)*patch.dx() + x_PmO;
  real y_p = (real(j) + 0.5)*patch.dy() + y_PmO;
  real z_p = (real(k) + 0.5)*patch.dz() + z_PmO;


  real p_sum = 0, w_sum = 0;
  bool updated = false;

  // if point l_im in cylinder, interpolate to l_p
  if (l_im > (3.5 - layer) && l_p < 0.5) {
    real weight;
    real dp;
    computeDpAndWeight(patch, x_PmO, y_PmO, z_PmO, x_p, y_p, z_p , i - 1, j, k, patch.dx(), dp, weight);
    p_sum += dp;
    w_sum += weight;
    updated = true;
  }
  // if point l_ip in cylinder, interpolate to l_p
  if (l_ip > (3.5 - layer) && l_p < 0.5) {
    real weight;
    real dp;
    computeDpAndWeight(patch, x_PmO, y_PmO, z_PmO, x_p, y_p, z_p , i + 1, j, k, patch.dx(), dp, weight);
    p_sum += dp;
    w_sum += weight;
    updated = true;
  }
  // if point l_jm in cylinder, interpolate to l_p
  if (l_jm > (3.5 - layer) && l_p < 0.5) {
    real weight;
    real dp;
    computeDpAndWeight(patch, x_PmO, y_PmO, z_PmO, x_p, y_p, z_p , i, j - 1, k, patch.dy(), dp, weight);
    p_sum += dp;
    w_sum += weight;
    updated = true;
  }
  // if point l_jp in cylinder, interpolate to l_p
  if (l_jp > (3.5 - layer) && l_p < 0.5) {
    real weight;
    real dp;
    computeDpAndWeight(patch, x_PmO, y_PmO, z_PmO, x_p, y_p, z_p , i, j + 1, k, patch.dy(), dp, weight);
    p_sum += dp;
    w_sum += weight;
    updated = true;
  }
  // if point l_km in cylinder, set zero gradient p
  if (l_km > (3.5 - layer) && l_p < 0.5) {
    real weight;
    real dp;
    computeDpAndWeight(patch, x_PmO, y_PmO, z_PmO, x_p, y_p, z_p , i, j, k - 1, 0, dp, weight);
    p_sum += dp;
    w_sum += weight;
    updated = true;
  }
  // if point l_kp in cylinder, set zero gradient p
  if (l_kp > (3.5 - layer) && l_p < 0.5) {
    real weight;
    real dp;
    computeDpAndWeight(patch, x_PmO, y_PmO, z_PmO, x_p, y_p, z_p , i, j, k + 1, 0, dp, weight);
    p_sum += dp;
    w_sum += weight;
    updated = true;
  }

  real var1[5];
  if (updated == true) {
    patch.f(2, 0, i, j, k) = 3 - layer;

    real u_layer =  omega*y_p;         // ==  m_Omega*norm(r_2)*r_2[1]/norm(r_2) ->  v_t*r*sin(th);
    real v_layer = -omega*x_p;         // == -m_Omega*norm(r_2)*r_2[0]/norm(r_2) -> -v_t*r*cos(th);

    PerfectGas::primitiveToConservative(p_sum/w_sum, t_inf, u_layer, v_layer, 0, var1);
    patch.setVar(dim, 0, i, j, k, var1);
    if ( 3 - layer > 1.5) {
    }
    else {
    }
  }
  else if (patch.f(2, 0, i, j, k) < 0.5) {
    PerfectGas::primitiveToConservative(p_inf, t_inf, 0, 0, 0, var1);
    patch.setVar(dim, 0, i, j, k, var1);
  }
}

void GPU_CylinderInCartisianPatch::operator ()()
{
  vec3_t x_PmO = m_Po - m_Xo;

  //size_t k_lines = max(size_t(1), size_t(m_MaxNumThreads/m_Patch.sizeK()));
  {
    //dim3 blocks(m_Patch.sizeI()/2+1, m_Patch.sizeJ()/k_lines+1, 1);
    //dim3 threads(m_Patch.sizeK(), k_lines, 1);

    dim3 blocks(m_Patch.sizeI()+1, m_Patch.sizeJ()+1, 1);
    dim3 threads(m_Patch.sizeK()+1, 1, 1);

    GPU_CylinderInCartisianPatch_Mark<<<blocks, threads>>>(m_Patch, x_PmO[0], x_PmO[1], x_PmO[2], m_Rad*m_Rad, m_Height*m_Height);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR;
    GPU_CylinderInCartisianPatch_ComputeLayer<<<blocks, threads>>>(m_Patch, x_PmO[0], x_PmO[1], x_PmO[2], m_Pinf, m_Tinf, m_Omega, 1);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR;
    GPU_CylinderInCartisianPatch_ComputeLayer<<<blocks, threads>>>(m_Patch, x_PmO[0], x_PmO[1], x_PmO[2], m_Pinf, m_Tinf, m_Omega, 2);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR;
  }

}

#endif // CUDA

#endif // GPU_CYLINDERINCARTISIANPATCH_H
