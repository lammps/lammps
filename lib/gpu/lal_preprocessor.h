// **************************************************************************
//                               preprocessor.h
//                             -------------------
//                           W. Michael Brown (ORNL)
//                           Nitin Dhamankar (Intel)
//
//  Device-side preprocessor definitions
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
//
//  For OpenCL, the configuration is a string (optionally controlled at
//  runtime) where tokens specify the values below in order)
//
//  CONFIG_ID:
//     Definition:   Unique ID for a configuration
//                   100-199 for NVIDIA GPUs with CUDA / HIP
//                   200-299 for NVIDIA GPUs with OpenCL
//                   300-399 for AMD GPUs with HIP
//                   400-499 for AMD GPUs with OpenCL
//                   500-599 for Intel GPUs with OpenCL
//  SIMD_SIZE:
//     Definition:   For CUDA this is the warp size.
//                   For AMD this is the wavefront size.
//                   For OpenCL < 2.1 this is the number of workitems
//                     guarenteed to have the same instruction pointer
//                   For OpenCL >= 2.1 this is the smallest expected subgroup
//                     size. Actually subgroup sizes are determined per kernel.
//  MEM_THREADS
//     Definition:   Number of elements in main memory transaction. Used in
//                   PPPM. If unknown, set to SIMD_SIZE.
//  SHUFFLE_AVAIL
//     Definition:   Controls the use of instructions for horizontal vector
//                   operations. 0 disables and will increase shared memory
//                   usage. 1 enables for CUDA, HIP, and OpenCL >= 2.1 on
//                   NVIDIA and Intel devices.
//  FAST_MATH
//     Definition:   0: do not use -cl-fast-relaxed-math optimization flag or
//                   native transcendentals for OpenCL (fused multiply-add
//                   still enabled). For CUDA and HIP, this is controlled by
//                   the Makefile at compile time. 1: enable fast math opts
//
//  THREADS_PER_ATOM
//     Definition:   Default number of work items or CUDA threads assigned per
//                   per atom for pair styles
//     Restrictions: Must be power of 2; THREADS_PER_ATOM<=SIMD_SIZE
//  THREADS_PER_CHARGE
//     Definition:   Default number of work items or CUDA threads assigned per
//                   per atom for pair styles using charge
//     Restrictions: Must be power of 2; THREADS_PER_ATOM<=SIMD_SIZE
//  THREADS_PER_THREE
//     Definition:   Default number of work items or CUDA threads assigned per
//                   per atom for 3-body styles
//     Restrictions: Must be power of 2; THREADS_PER_ATOM^2<=SIMD_SIZE
//
//  BLOCK_PAIR
//     Definition:   Default block size for pair styles
//     Restrictions: Must be integer multiple of SIMD_SIZE
//  BLOCK_BIO_PAIR
//     Definition:   Default block size for CHARMM styles
//     Restrictions: Must be integer multiple of SIMD_SIZE
//  BLOCK_ELLIPSE
//     Definition:   Default block size for ellipsoidal models and some 3-body
//                   styles
//     Restrictions: Must be integer multiple of SIMD_SIZE
//  PPPM_BLOCK_1D
//     Definition:   Default block size for PPPM kernels
//     Restrictions: Must be integer multiple of SIMD_SIZE
//  BLOCK_NBOR_BUILD
//     Definition:   Default block size for neighbor list builds
//     Restrictions: Must be integer multiple of SIMD_SIZE
//  BLOCK_CELL_2D
//     Definition:   Default block size in each dimension for matrix transpose
//  BLOCK_CELL_ID
//     Definition:   Unused in current implementation; Maintained for legacy
//                   purposes and specialized builds
//
//  MAX_SHARED_TYPES 8
//     Definition:   Max # of atom type params can be stored in shared memory
//     Restrictions: MAX_SHARED_TYPES*MAX_SHARED_TYPES<=BLOCK_PAIR
//  MAX_BIO_SHARED_TYPES
//     Definition:   Max # of atom type params can be stored in shared memory
//     Restrictions: MAX_BIO_SHARED_TYPES<=BLOCK_BIO_PAIR*2
//  PPPM_MAX_SPLINE
//     Definition:   Maximum order for splines in PPPM
//     Restrictions: PPPM_BLOCK_1D>=PPPM_MAX_SPLINE*PPPM_MAX_SPLINE
//
//*************************************************************************/

// -------------------------------------------------------------------------
//                           CUDA and HIP DEFINITIONS
// -------------------------------------------------------------------------

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_pre_cuda_hip.h"
#define ucl_pow pow
#endif

// -------------------------------------------------------------------------
//                         OPENCL DEVICE CONFIGURATAIONS
// -------------------------------------------------------------------------

// See lal_pre_ocl_config.h for OpenCL device configurations

#if !defined(NV_KERNEL) && !defined(USE_HIP)

#define USE_OPENCL

// -------------------------------------------------------------------------
//                         OPENCL KERNEL MACROS
// -------------------------------------------------------------------------

#if (__OPENCL_VERSION__ > 199)
#define NOUNROLL __attribute__((opencl_unroll_hint(1)))
#else
#define NOUNROLL
#endif

#define GLOBAL_ID_X get_global_id(0)
#define THREAD_ID_X get_local_id(0)
#define BLOCK_ID_X get_group_id(0)
#define BLOCK_SIZE_X get_local_size(0)
#define GLOBAL_SIZE_X get_global_size(0)
#define THREAD_ID_Y get_local_id(1)
#define BLOCK_ID_Y get_group_id(1)
#define NUM_BLOCKS_X get_num_groups(0)
#define __syncthreads() barrier(CLK_LOCAL_MEM_FENCE)
#define ucl_inline inline

// -------------------------------------------------------------------------
//                      OPENCL KERNEL MACROS - TEXTURES
// -------------------------------------------------------------------------

#define fetch4(ans,i,x) ans=x[i]
#define fetch(ans,i,q) ans=q[i]

// -------------------------------------------------------------------------
//                       OPENCL KERNEL MACROS - MATH
// -------------------------------------------------------------------------

#ifndef _SINGLE_SINGLE

#ifndef cl_khr_fp64
#ifndef cl_amd_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
#endif
#if defined(cl_khr_fp64)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif

#endif

#define fast_mul(X,Y) (X)*(Y)

#define ucl_atan atan
#define ucl_cbrt cbrt
#define ucl_ceil ceil
#define ucl_abs fabs
#define ucl_erfc erfc

#if defined(FAST_MATH) && !defined(_DOUBLE_DOUBLE)

#define ucl_exp native_exp
#define ucl_pow pow
#define ucl_powr native_powr
#define ucl_rsqrt native_rsqrt
#define ucl_sqrt native_sqrt
#define ucl_recip native_recip

#else

#define ucl_exp exp
#define ucl_pow pow
#define ucl_powr powr
#define ucl_rsqrt rsqrt
#define ucl_sqrt sqrt
#define ucl_recip(x) ((numtyp)1.0/(x))

#endif

// -------------------------------------------------------------------------
//                      OPENCL KERNEL MACROS - SHUFFLE
// -------------------------------------------------------------------------

#if (SHUFFLE_AVAIL == 1)
  #ifdef cl_intel_subgroups
    #pragma OPENCL EXTENSION cl_intel_subgroups : enable
    #define shfl_down(var, delta, width) \
      intel_sub_group_shuffle_down(var, var, delta)
    #define shfl_xor(var, lanemask, width) \
      intel_sub_group_shuffle_xor(var, lanemask)
    #define simd_broadcast_i(var, src, width) sub_group_broadcast(var, src)
    #define simd_broadcast_f(var, src, width) sub_group_broadcast(var, src)
    #define simd_broadcast_d(var, src, width) sub_group_broadcast(var, src)
  #else
    #ifdef _SINGLE_SINGLE
      inline float shfl_down(float var, unsigned int delta, int width) {
        float ret;
        int c;
        c = ((SIMD_SIZE-width) << 8) | 0x1f;
        asm volatile ("shfl.sync.down.b32 %0, %1, %2, %3, %4;" : "=f"(ret) : "f"(var), "r"(delta), "r"(c), "r"(0xffffffff));
        return ret;
      }
      inline float shfl_xor(float var, unsigned int lanemask, int width) {
        float ret;
        int c;
        c = ((SIMD_SIZE-width) << 8) | 0x1f;
        asm volatile ("shfl.sync.bfly.b32 %0, %1, %2, %3, %4;" : "=f"(ret) : "f"(var), "r"(lanemask), "r"(c), "r"(0xffffffff));
        return ret;
      }
    #else
      inline double shfl_down(double var, unsigned int delta, int width) {
        int c = ((SIMD_SIZE-width) << 8) | 0x1f;
        int x,y,x2,y2;
        double ans;
        asm volatile ("mov.b64 {%0, %1}, %2;" : "=r"(y), "=r"(x) : "d"(var));
        asm volatile ("shfl.sync.down.b32 %0, %1, %2, %3, %4;" : "=r"(x2) : "r"(x), "r"(delta), "r"(c), "r"(0xffffffff));
        asm volatile ("shfl.sync.down.b32 %0, %1, %2, %3, %4;" : "=r"(y2) : "r"(y), "r"(delta), "r"(c), "r"(0xffffffff));
        asm volatile ("mov.b64 %0, {%1, %2};" : "=d"(ans) : "r"(y2), "r"(x2));
        return ans;
      }
      inline double shfl_xor(double var, unsigned int lanemask, int width) {
        int c = ((SIMD_SIZE-width) << 8) | 0x1f;
        int x,y,x2,y2;
        double ans;
        asm volatile ("mov.b64 {%0, %1}, %2;" : "=r"(y), "=r"(x) : "d"(var));
        asm volatile ("shfl.sync.bfly.b32 %0, %1, %2, %3, %4;" : "=r"(x2) : "r"(x), "r"(lanemask), "r"(c), "r"(0xffffffff));
        asm volatile ("shfl.sync.bfly.b32 %0, %1, %2, %3, %4;" : "=r"(y2) : "r"(y), "r"(lanemask), "r"(c), "r"(0xffffffff));
        asm volatile ("mov.b64 %0, {%1, %2};" : "=d"(ans) : "r"(y2), "r"(x2));
        return ans;
      }
    #endif
    inline int simd_broadcast_i(int var, unsigned int src, int width) {
      int ret;
      int c;
      c = ((SIMD_SIZE-width) << 8) | 0x1f;
      asm volatile ("shfl.sync.idx.b32 %0, %1, %2, %3, %4;" : "=f"(ret) : "f"(var), "r"(src), "r"(c), "r"(0xffffffff));
      return ret;
    }
    inline float simd_broadcast_f(float var, unsigned int src, int width) {
      float ret;
      int c;
      c = ((SIMD_SIZE-width) << 8) | 0x1f;
      asm volatile ("shfl.sync.idx.b32 %0, %1, %2, %3, %4;" : "=f"(ret) : "f"(var), "r"(src), "r"(c), "r"(0xffffffff));
      return ret;
    }
    #ifdef _DOUBLE_DOUBLE
      inline double simd_broadcast_d(double var, unsigned int src, int width) {
        int c = ((SIMD_SIZE-width) << 8) | 0x1f;
        int x,y,x2,y2;
        double ans;
        asm volatile ("mov.b64 {%0, %1}, %2;" : "=r"(y), "=r"(x) : "d"(var));
        asm volatile ("shfl.sync.idx.b32 %0, %1, %2, %3, %4;" : "=r"(x2) : "r"(x), "r"(src), "r"(c), "r"(0xffffffff));
        asm volatile ("shfl.sync.idx.b32 %0, %1, %2, %3, %4;" : "=r"(y2) : "r"(y), "r"(src), "r"(c), "r"(0xffffffff));
        asm volatile ("mov.b64 %0, {%1, %2};" : "=d"(ans) : "r"(y2), "r"(x2));
        return ans;
      }
    #endif
  #endif
#endif

// -------------------------------------------------------------------------
//                      OPENCL KERNEL MACROS - SUBGROUPS
// -------------------------------------------------------------------------

#ifdef USE_OPENCL_SUBGROUPS
  #ifndef cl_intel_subgroups
    #pragma OPENCL EXTENSION cl_khr_subgroups : enable
  #endif
  #define simdsync() sub_group_barrier(CLK_LOCAL_MEM_FENCE)
  #define simd_size() get_max_sub_group_size()
#else
  #define simdsync()
  #define simd_size() SIMD_SIZE
#endif

// -------------------------------------------------------------------------
//                            END OPENCL DEFINITIONS
// -------------------------------------------------------------------------

#endif

// -------------------------------------------------------------------------
//                  ARCHITECTURE INDEPENDENT DEFINITIONS
// -------------------------------------------------------------------------

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp2 double2
#define numtyp4 double4
#define acctyp double
#define acctyp2 double2
#define acctyp4 double4
#endif

#ifdef _SINGLE_DOUBLE
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp double
#define acctyp2 double2
#define acctyp4 double4
#endif

#ifndef numtyp
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp float
#define acctyp2 float2
#define acctyp4 float4
#endif

#define EWALD_F (numtyp)1.12837917
#define EWALD_P (numtyp)0.3275911
#define A1 (numtyp)0.254829592
#define A2 (numtyp)-0.284496736
#define A3 (numtyp)1.421413741
#define A4 (numtyp)-1.453152027
#define A5 (numtyp)1.061405429

#define SBBITS 30
#define NEIGHMASK 0x3FFFFFFF
ucl_inline int sbmask(int j) { return j >> SBBITS & 3; };

#define SBBITS15 29
#define NEIGHMASK15 0x1FFFFFFF
ucl_inline int sbmask15(int j) { return j >> SBBITS15 & 7; };

// default to 32-bit smallint and other ints, 64-bit bigint:
// same as defined in src/lmptype.h
#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) && \
    !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif
