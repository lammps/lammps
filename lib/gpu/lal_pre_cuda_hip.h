// **************************************************************************
//                               pre_cuda_hip.h
//                             -------------------
//                           W. Michael Brown (ORNL)
//                           Nitin Dhamankar (Intel)
//
//  Device-side preprocessor definitions for CUDA and HIP builds
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************/

//*************************************************************************
//                       Device Configuration Definitions
//                    See lal_preprocessor.h for definitions
//*************************************************************************/

// -------------------------------------------------------------------------
//                           CUDA and HIP DEFINITIONS
// -------------------------------------------------------------------------

#if defined(NV_KERNEL) || defined(USE_HIP)

// -------------------------------------------------------------------------
//                             DEVICE CONFIGURATION
// -------------------------------------------------------------------------


#if defined(__HIP_PLATFORM_HCC__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIP_PLATFORM_SPIRV__)
#define CONFIG_ID 303
#define SIMD_SIZE 64
#else
#define CONFIG_ID 103
#define SIMD_SIZE 32
#endif

#define MEM_THREADS SIMD_SIZE
#define SHUFFLE_AVAIL 1
#define FAST_MATH 1

#define THREADS_PER_ATOM 4
#define THREADS_PER_CHARGE 8
#define THREADS_PER_THREE 2

#define BLOCK_PAIR 256
#define BLOCK_BIO_PAIR 256
#define BLOCK_ELLIPSE 128
#define PPPM_BLOCK_1D 64
#define BLOCK_NBOR_BUILD 128
#define BLOCK_CELL_2D 8
#define BLOCK_CELL_ID 128

#define MAX_SHARED_TYPES 11
#define MAX_BIO_SHARED_TYPES 128
#define PPPM_MAX_SPLINE 8

// -------------------------------------------------------------------------
//                              KERNEL MACROS
// -------------------------------------------------------------------------

#ifdef USE_HIP
#include <hip/hip_runtime.h>
#endif

#define fast_mul(X,Y) (X)*(Y)

#define EVFLAG 1
#define NOUNROLL
#define GLOBAL_ID_X threadIdx.x+fast_mul(blockIdx.x,blockDim.x)
#define GLOBAL_ID_Y threadIdx.y+fast_mul(blockIdx.y,blockDim.y)
#define GLOBAL_SIZE_X fast_mul(gridDim.x,blockDim.x);
#define GLOBAL_SIZE_Y fast_mul(gridDim.y,blockDim.y);
#define THREAD_ID_X threadIdx.x
#define THREAD_ID_Y threadIdx.y
#define BLOCK_ID_X blockIdx.x
#define BLOCK_ID_Y blockIdx.y
#define BLOCK_SIZE_X blockDim.x
#define BLOCK_SIZE_Y blockDim.y
#define NUM_BLOCKS_X gridDim.x

#define __kernel extern "C" __global__
#ifdef __local
#undef __local
#endif
#define __local __shared__
#define __global
#define restrict __restrict__
#define atom_add atomicAdd
#define ucl_inline static __inline__ __device__

#define simd_size() SIMD_SIZE

#define simdsync()

#ifdef NV_KERNEL
#if (__CUDACC_VER_MAJOR__ >= 9)
#undef simdsync
#define simdsync() __syncwarp(0xffffffff)
#endif
#endif

#ifdef __HIP_PLATFORM_NVCC__
#undef simdsync()
#define simdsync() __syncwarp(0xffffffff)
#endif

// -------------------------------------------------------------------------
//                         KERNEL MACROS - TEXTURES
// -------------------------------------------------------------------------

#if defined(__HIP_PLATFORM_HCC__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIP_PLATFORM_SPIRV__)
#define _texture(name, type)  __device__ type* name
#define _texture_2d(name, type)  __device__ type* name
#else
#define _texture(name, type)  texture<type> name
#define _texture_2d(name, type) texture<type,1> name
#endif

#if (__CUDACC_VER_MAJOR__ < 11)
  #ifdef _DOUBLE_DOUBLE
  #define fetch4(ans,i,pos_tex) {                        \
    int4 xy = tex1Dfetch(pos_tex,i*2);                   \
    int4 zt = tex1Dfetch(pos_tex,i*2+1);                 \
    ans.x=__hiloint2double(xy.y, xy.x);                  \
    ans.y=__hiloint2double(xy.w, xy.z);                  \
    ans.z=__hiloint2double(zt.y, zt.x);                  \
    ans.w=__hiloint2double(zt.w, zt.z);                  \
  }
  #define fetch(ans,i,q_tex) {                           \
    int2 qt = tex1Dfetch(q_tex,i);                       \
    ans=__hiloint2double(qt.y, qt.x);                    \
  }
  #elif  defined(__HIP_PLATFORM_SPIRV__)
      #define fetch4(ans,i,pos_tex) tex1Dfetch(&ans, pos_tex, i);
      #define fetch(ans,i,q_tex) tex1Dfetch(&ans, q_tex,i);
  #else
    #define fetch4(ans,i,pos_tex) ans=tex1Dfetch(pos_tex, i);
    #define fetch(ans,i,q_tex) ans=tex1Dfetch(q_tex,i);
  #endif
#else
  #define fetch4(ans,i,x) ans=x[i]
  #define fetch(ans,i,q) ans=q[i]
  #undef _texture
  #undef _texture_2d
  #define _texture(name, type)
  #define _texture_2d(name, type)
  #define pos_tex x_
  #define quat_tex qif
  #define q_tex q_
  #define vel_tex v_
  #define mu_tex mu_
#endif

#if defined(__HIP_PLATFORM_HCC__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIP_PLATFORM_SPIRV__)

#undef fetch4
#undef fetch

#ifdef _DOUBLE_DOUBLE
#define fetch4(ans,i,pos_tex) (ans=*(((double4*)pos_tex) + i))
#define fetch(ans,i,q_tex)    (ans=*(((double *)  q_tex) + i))
#else
#define fetch4(ans,i,pos_tex) (ans=*(((float4*)pos_tex) + i))
#define fetch(ans,i,q_tex)    (ans=*(((float *)  q_tex) + i))
#endif

#endif

// -------------------------------------------------------------------------
//                           KERNEL MACROS - MATH
// -------------------------------------------------------------------------

#ifdef _DOUBLE_DOUBLE

#define ucl_exp exp
#define ucl_powr pow
#define ucl_atan atan
#define ucl_cbrt cbrt
#define ucl_ceil ceil
#define ucl_abs fabs
#define ucl_recip(x) ((numtyp)1.0/(x))
#define ucl_rsqrt rsqrt
#define ucl_sqrt sqrt
#define ucl_erfc erfc

#else

#define ucl_exp expf
#define ucl_powr powf
#define ucl_atan atanf
#define ucl_cbrt cbrtf
#define ucl_ceil ceilf
#define ucl_abs fabsf
#define ucl_recip(x) ((numtyp)1.0/(x))
#define ucl_rsqrt rsqrtf
#define ucl_sqrt sqrtf
#define ucl_erfc erfcf

#endif

// -------------------------------------------------------------------------
//                         KERNEL MACROS - SHUFFLE
// -------------------------------------------------------------------------

#if SHUFFLE_AVAIL == 1

#ifndef USE_HIP
#if (__CUDACC_VER_MAJOR__ < 9)
#define CUDA_PRE_NINE
#endif
#endif

#if defined(CUDA_PRE_NINE) || defined(__HIP_PLATFORM_HCC__) || defined(__HIP_PLATFORM_AMD__) || defined(__HIP_PLATFORM_SPIRV__)

  #ifdef _SINGLE_SINGLE
    #define shfl_down __shfl_down
    #define shfl_xor __shfl_xor
  #else
    ucl_inline double shfl_down(double var, unsigned int delta, int width) {
      int2 tmp;
      tmp.x = __double2hiint(var);
      tmp.y = __double2loint(var);
      tmp.x = __shfl_down(tmp.x,delta,width);
      tmp.y = __shfl_down(tmp.y,delta,width);
      return __hiloint2double(tmp.x,tmp.y);
    }
    ucl_inline double shfl_xor(double var, unsigned int lanemask, int width) {
      int2 tmp;
      tmp.x = __double2hiint(var);
      tmp.y = __double2loint(var);
      tmp.x = __shfl_xor(tmp.x,lanemask,width);
      tmp.y = __shfl_xor(tmp.y,lanemask,width);
      return __hiloint2double(tmp.x,tmp.y);
    }
  #endif
  #define simd_broadcast_i __shfl
  #define simd_broadcast_f __shfl
  #ifdef _DOUBLE_DOUBLE
    ucl_inline double simd_broadcast_d(double var, unsigned int src,
                                       int width) {
      int2 tmp;
      tmp.x = __double2hiint(var);
      tmp.y = __double2loint(var);
      tmp.x = __shfl(tmp.x,src,width);
      tmp.y = __shfl(tmp.y,src,width);
      return __hiloint2double(tmp.x,tmp.y);
    }
  #endif

#else

  #ifdef _SINGLE_SINGLE
  ucl_inline float shfl_down(float var, unsigned int delta, int width) {
    return __shfl_down_sync(0xffffffff, var, delta, width);
  }
  ucl_inline float shfl_xor(float var, unsigned int lanemask, int width) {
    return __shfl_xor_sync(0xffffffff, var, lanemask, width);
  }
  #else
  ucl_inline double shfl_down(double var, unsigned int delta, int width) {
    int2 tmp;
    tmp.x = __double2hiint(var);
    tmp.y = __double2loint(var);
    tmp.x = __shfl_down_sync(0xffffffff,tmp.x,delta,width);
    tmp.y = __shfl_down_sync(0xffffffff,tmp.y,delta,width);
    return __hiloint2double(tmp.x,tmp.y);
  }
  ucl_inline double shfl_xor(double var, unsigned int lanemask, int width) {
    int2 tmp;
    tmp.x = __double2hiint(var);
    tmp.y = __double2loint(var);
    tmp.x = __shfl_xor_sync(0xffffffff,tmp.x,lanemask,width);
    tmp.y = __shfl_xor_sync(0xffffffff,tmp.y,lanemask,width);
    return __hiloint2double(tmp.x,tmp.y);
  }
  #endif
  #define simd_broadcast_i(var, src, width) \
    __shfl_sync(0xffffffff, var, src, width)
  #define simd_broadcast_f(var, src, width) \
    __shfl_sync(0xffffffff, var, src, width)
  #ifdef _DOUBLE_DOUBLE
  ucl_inline double simd_broadcast_d(double var, unsigned int src, int width) {
    int2 tmp;
    tmp.x = __double2hiint(var);
    tmp.y = __double2loint(var);
    tmp.x = __shfl_sync(0xffffffff,tmp.x,src,width);
    tmp.y = __shfl_sync(0xffffffff,tmp.y,src,width);
    return __hiloint2double(tmp.x,tmp.y);
  }
  #endif
#endif

#endif

// -------------------------------------------------------------------------
//                            END CUDA / HIP DEFINITIONS
// -------------------------------------------------------------------------

#endif
