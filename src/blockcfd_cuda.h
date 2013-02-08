#ifndef BLOCKCFD_CUDA_H
#define BLOCKCFD_CUDA_H

#ifdef __CUDACC__
  #define CUDA_DO __device__
  #define CUDA_HO __host__
  #define CUDA_DH __device__ __host__
#else
  #define CUDA_DO
  #define CUDA_HO
  #define CUDA_DH
#endif

#endif // BLOCKCFD_CUDA_H
