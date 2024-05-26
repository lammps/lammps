// clang-format off
/* -*- c++ -*- -------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "lmptype.h"

#ifdef __INTEL_LLVM_COMPILER
#define USE_OMP_SIMD
#define __INTEL_COMPILER __INTEL_LLVM_COMPILER
#define __INTEL_COMPILER_BUILD_DATE __INTEL_LLVM_COMPILER
// Indicate to vectorizer that it is safe to use dword indexed gather
#define IP_PRE_dword_index(i) ((i) & NEIGHMASK)
#else
#define IP_PRE_dword_index(i) i
#endif

#ifdef __INTEL_COMPILER
#define LMP_SIMD_COMPILER
#if (__INTEL_COMPILER_BUILD_DATE > 20160720)
#define LMP_INTEL_USE_SIMDOFF
#endif
#pragma warning (disable:3948)
#pragma warning (disable:3949)
#pragma warning (disable:13200)
#endif

#ifdef __INTEL_OFFLOAD
#ifdef LMP_INTEL_OFFLOAD
#define _LMP_INTEL_OFFLOAD
#ifdef __TARGET_ARCH_MIC
#ifndef __MIC__
#define __MIC__ 1
#endif
#endif
#endif
#endif

#ifndef LMP_INTEL_PREPROCESS_H
#define LMP_INTEL_PREPROCESS_H

// LAMMPS_MEMALIGN is set to 64 by default for -DLMP_INTEL
// so we only need to error out in case of a different alignment
#if LAMMPS_MEMALIGN && (LAMMPS_MEMALIGN != 64)
#error Please set -DLAMMPS_MEMALIGN=64 in CCFLAGS of your LAMMPS makefile for INTEL package
#endif

#if defined(_OPENMP)
#define _use_omp_pragma(txt) _Pragma(txt)
#else
#define _use_omp_pragma(txt)
#endif

#if defined(LMP_SIMD_COMPILER)
#define _use_simd_pragma(txt) _Pragma(txt)
#else
#define _use_simd_pragma(txt)
#endif

namespace LAMMPS_NS {

enum {LMP_OVERFLOW, LMP_LOCAL_MIN, LMP_LOCAL_MAX, LMP_GHOST_MIN,
      LMP_GHOST_MAX};
enum {TIME_PACK, TIME_HOST_NEIGHBOR, TIME_HOST_PAIR, TIME_OFFLOAD_NEIGHBOR,
      TIME_OFFLOAD_PAIR, TIME_OFFLOAD_WAIT, TIME_OFFLOAD_LATENCY,
      TIME_IMBALANCE};

#define NUM_ITIMERS ( TIME_IMBALANCE + 1 )
#define INTEL_MIC_VECTOR_WIDTH 16
#define INTEL_VECTOR_WIDTH 4
#define INTEL_MAX_STENCIL 256
// INTEL_MAX_STENCIL * sqrt(INTEL_MAX_STENCIL)
#define INTEL_MAX_STENCIL_CHECK 4096
#define INTEL_P3M_MAXORDER 8
#define INTEL_P3M_ALIGNED_MAXORDER 8

#ifdef __INTEL_COMPILER
#ifdef __AVX__
#undef INTEL_VECTOR_WIDTH
#define INTEL_VECTOR_WIDTH 8
#endif

#ifdef __AVX2__
#undef INTEL_VECTOR_WIDTH
#define INTEL_VECTOR_WIDTH 8
#endif

#ifdef __AVX512F__
#undef INTEL_VECTOR_WIDTH
#define INTEL_VECTOR_WIDTH 16
#define INTEL_V512 1
#define INTEL_VMASK 1
#else
#ifdef __MIC__
#define INTEL_V512 1
#define INTEL_VMASK 1
#define INTEL_HTHREADS 4
#endif
#endif

#ifdef __AVX512ER__
#define INTEL_HTHREADS 4
#endif

#ifdef __AVX512CD__
#ifndef _LMP_INTEL_OFFLOAD
#define LMP_USE_AVXCD
#endif
#endif

#ifdef __MIC__
#define INTEL_COMPILE_WIDTH INTEL_MIC_VECTOR_WIDTH
#else
#define INTEL_COMPILE_WIDTH INTEL_VECTOR_WIDTH
#endif

#else

#undef INTEL_VECTOR_WIDTH
#define INTEL_VECTOR_WIDTH 1
#define INTEL_COMPILE_WIDTH 1
#if defined(__AVX512F__) && !defined(INTEL_VMASK)
#define INTEL_VMASK 1
#endif

#endif

#define INTEL_DATA_ALIGN 64
#define INTEL_ONEATOM_FACTOR 1
#define INTEL_MIC_NBOR_PAD INTEL_MIC_VECTOR_WIDTH
#define INTEL_NBOR_PAD INTEL_VECTOR_WIDTH
#define INTEL_LB_MEAN_WEIGHT 0.1
#define INTEL_BIGP 1e15
#define INTEL_MAX_HOST_CORE_COUNT 512
#define INTEL_MAX_COI_CORES 36

#ifndef INTEL_HTHREADS
#define INTEL_HTHREADS 2
#endif

#if INTEL_DATA_ALIGN > 1

#define IP_PRE_edge_align(n, esize)                                     \
  {                                                                     \
    const int pad_mask = ~static_cast<int>(INTEL_DATA_ALIGN/esize-1);   \
    n = (n + INTEL_DATA_ALIGN / esize - 1) & pad_mask;                  \
  }

#else

#define IP_PRE_edge_align(n, esize)                                     \

#endif

#define IP_PRE_get_stride(stride, n, datasize, torque)          \
  {                                                             \
    int blength = n;                                            \
    if (torque) blength *= 2;                                   \
    const int bytes = blength * datasize;                       \
    stride = INTEL_DATA_ALIGN - (bytes % INTEL_DATA_ALIGN);     \
    stride = blength + stride / datasize;                       \
  }

#if defined(_OPENMP)

#define IP_PRE_omp_range(ifrom, ito, tid, inum, nthreads)       \
  {                                                             \
    int idelta = inum/nthreads;                                 \
    const int imod = inum % nthreads;                           \
    ifrom = tid * idelta;                                       \
    ito = ifrom + idelta;                                       \
    if (tid < imod) {                                           \
      ito+=tid+1;                                               \
      ifrom+=tid;                                               \
    } else {                                                    \
      ito+=imod;                                                \
      ifrom+=imod;                                              \
    }                                                           \
  }

#define IP_PRE_omp_range_id(ifrom, ito, tid, inum, nthreads)    \
  {                                                             \
    tid = omp_get_thread_num();                                 \
    IP_PRE_omp_range(ifrom, ito, tid, inum, nthreads);          \
  }

#define IP_PRE_omp_stride(ifrom, ip, ito, tid, inum, nthr)      \
  {                                                             \
    if (nthr <= INTEL_HTHREADS) {                               \
      ifrom = tid;                                              \
      ito = inum;                                               \
      ip = nthr;                                                \
    } else if (nthr % INTEL_HTHREADS == 0) {                    \
      int nd = nthr / INTEL_HTHREADS;                           \
      int td = tid / INTEL_HTHREADS;                            \
      int tm = tid % INTEL_HTHREADS;                            \
      IP_PRE_omp_range(ifrom, ito, td, inum, nd);               \
      ifrom += tm;                                              \
      ip = INTEL_HTHREADS;                                      \
    } else {                                                    \
      IP_PRE_omp_range(ifrom, ito, tid, inum, nthr);            \
      ip = 1;                                                   \
    }                                                           \
  }

#define IP_PRE_omp_stride_id(ifrom, ip, ito, tid, inum, nthr)   \
  {                                                             \
    tid = omp_get_thread_num();                                 \
    IP_PRE_omp_stride(ifrom, ip, ito, tid, inum, nthr);         \
  }

#define IP_PRE_omp_range_align(ifrom, ito, tid, inum, nthreads, \
                             datasize)                          \
{                                                               \
  int chunk_size = INTEL_DATA_ALIGN / datasize;                 \
  int idelta = static_cast<int>(ceil(static_cast<float>(inum)   \
                                     /chunk_size/nthreads));    \
  idelta *= chunk_size;                                         \
  ifrom = tid*idelta;                                           \
  ito = ifrom + idelta;                                         \
  if (ito > inum) ito = inum;                                   \
}

#define IP_PRE_omp_range_id_align(ifrom, ito, tid, inum,        \
                                nthreads, datasize)             \
  {                                                             \
    tid = omp_get_thread_num();                                 \
    IP_PRE_omp_range_align(ifrom, ito, tid, inum, nthreads,     \
                           datasize);                           \
  }

#define IP_PRE_omp_range_vec(ifrom, ito, tid, inum, nthreads,   \
                             vecsize)                           \
  {                                                             \
    int idelta = static_cast<int>(ceil(static_cast<float>(inum) \
                                       /vecsize/nthreads));     \
    idelta *= vecsize;                                          \
    ifrom = tid*idelta;                                         \
    ito = ifrom + idelta;                                       \
    if (ito > inum) ito = inum;                                 \
  }

#define IP_PRE_omp_range_id_vec(ifrom, ito, tid, inum,          \
                                nthreads, vecsize)              \
  {                                                             \
    tid = omp_get_thread_num();                                 \
    IP_PRE_omp_range_vec(ifrom, ito, tid, inum, nthreads,       \
                         vecsize);                              \
  }

#define IP_PRE_omp_stride_id_vec(ifrom, ip, ito, tid, inum,     \
                                 nthr, vecsize)                 \
  {                                                             \
    tid = omp_get_thread_num();                                 \
    if (nthr <= INTEL_HTHREADS) {                               \
      ifrom = tid*vecsize;                                      \
      ito = inum;                                               \
      ip = nthr*vecsize;                                        \
    } else if (nthr % INTEL_HTHREADS == 0) {                    \
      int nd = nthr / INTEL_HTHREADS;                           \
      int td = tid / INTEL_HTHREADS;                            \
      int tm = tid % INTEL_HTHREADS;                            \
      IP_PRE_omp_range_vec(ifrom, ito, td, inum, nd, vecsize);  \
      ifrom += tm * vecsize;                                    \
      ip = INTEL_HTHREADS * vecsize;                            \
    } else {                                                    \
      IP_PRE_omp_range_vec(ifrom, ito, tid, inum, nthr,         \
                           vecsize);                            \
      ip = vecsize;                                             \
    }                                                           \
  }

#else

#define IP_PRE_omp_range_id(ifrom, ito, tid, inum, nthreads)    \
  {                                                             \
    tid = 0;                                                    \
    ifrom = 0;                                                  \
    ito = inum;                                                 \
  }

#define IP_PRE_omp_range(ifrom, ito, tid, inum, nthreads)       \
  {                                                             \
    ifrom = 0;                                                  \
    ito = inum;                                                 \
  }

#define IP_PRE_omp_stride_id(ifrom, ip, ito, tid, inum, nthr)   \
  {                                                             \
    tid = 0;                                                    \
    ifrom = 0;                                                  \
    ito = inum;                                                 \
    ip = 1;                                                     \
  }

#define IP_PRE_omp_range_align(ifrom, ito, tid, inum, nthreads, \
                             datasize)                          \
{                                                               \
    ifrom = 0;                                                  \
    ito = inum;                                                 \
}

#define IP_PRE_omp_range_id_align(ifrom, ito, tid, inum,        \
                                nthreads, datasize)             \
{                                                               \
  tid = 0;                                                      \
  ifrom = 0;                                                    \
  ito = inum;                                                   \
}

#define IP_PRE_omp_range_id_vec(ifrom, ito, tid, inum,          \
                                nthreads, vecsize)              \
  {                                                             \
    tid = 0;                                                    \
    ifrom = 0;                                                  \
    ito = inum;                                                 \
  }

#define IP_PRE_omp_stride_id_vec(ifrom, ip, ito, tid, inum,     \
                                 nthr, vecsize)                 \
  {                                                             \
    tid = 0;                                                    \
    ifrom = 0;                                                  \
    ip = vecsize;                                               \
    ito = inum;                                                 \
  }

#endif

// TO BE DEPRECATED
#ifndef USE_OMP_SIMD

#define IP_PRE_fdotr_acc_force_l5(lf, lt, minlocal, nthreads, f_start,  \
                                  f_stride, pos, ov0, ov1, ov2,         \
                                  ov3, ov4, ov5)                        \
{                                                                       \
  acc_t *f_scalar = &f_start[0].x;                                      \
  flt_t *x_scalar = &pos[minlocal].x;                                   \
  int f_stride4 = f_stride * 4;                                         \
  _alignvar(acc_t ovv[16],64);                                          \
  int vwidth;                                                           \
  if (sizeof(acc_t) == sizeof(double))                                  \
    vwidth = INTEL_COMPILE_WIDTH/2;                                     \
  else                                                                  \
    vwidth = INTEL_COMPILE_WIDTH;                                       \
  if (vwidth < 4) vwidth = 4;                                           \
  _use_simd_pragma("vector aligned")                                    \
  _use_simd_pragma("simd")                                              \
  for (int v = 0; v < vwidth; v++) ovv[v] = (acc_t)0.0;                 \
  int remainder = lt % vwidth;                                          \
  if (lf > lt) remainder = 0;                                           \
  const int v_range = lt - remainder;                                   \
  if (nthreads == 2) {                                                  \
    acc_t *f_scalar2 = f_scalar + f_stride4;                            \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("vector aligned")                                \
      _use_simd_pragma("simd")                                          \
      for (int v = 0; v < vwidth; v++) {                                \
        f_scalar[n+v] += f_scalar2[n+v];                                \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      }                                                                 \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    _use_simd_pragma("loop_count min(4) max(INTEL_COMPILE_WIDTH)")      \
    for (int n = v_range; n < lt; n++)                                  \
      f_scalar[n] += f_scalar2[n];                                      \
  } else if (nthreads==4) {                                             \
    acc_t *f_scalar2 = f_scalar + f_stride4;                            \
    acc_t *f_scalar3 = f_scalar2 + f_stride4;                           \
    acc_t *f_scalar4 = f_scalar3 + f_stride4;                           \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("vector aligned")                                \
      _use_simd_pragma("simd")                                          \
      for (int v = 0; v < vwidth; v++) {                                \
        f_scalar[n+v] += f_scalar2[n+v] + f_scalar3[n+v] +              \
          f_scalar4[n+v];                                               \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      }                                                                 \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    _use_simd_pragma("loop_count min(4) max(INTEL_COMPILE_WIDTH)")      \
    for (int n = v_range; n < lt; n++)                                  \
      f_scalar[n] += f_scalar2[n] + f_scalar3[n] + f_scalar4[n];        \
  } else if (nthreads==1) {                                             \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("vector aligned")                                \
      _use_simd_pragma("simd")                                          \
      for (int v = 0; v < vwidth; v++)                                  \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
  } else if (nthreads==3) {                                             \
    acc_t *f_scalar2 = f_scalar + f_stride4;                            \
    acc_t *f_scalar3 = f_scalar2 + f_stride4;                           \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("vector aligned")                                \
      _use_simd_pragma("simd")                                          \
      for (int v = 0; v < vwidth; v++) {                                \
        f_scalar[n+v] += f_scalar2[n+v] + f_scalar3[n+v];               \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      }                                                                 \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    _use_simd_pragma("loop_count min(4) max(INTEL_COMPILE_WIDTH)")      \
    for (int n = v_range; n < lt; n++)                                  \
      f_scalar[n] += f_scalar2[n] + f_scalar3[n];                       \
  }                                                                     \
  for (int n = v_range; n < lt; n += 4) {                               \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    for (int v = 0; v < 4; v++)                                         \
      ovv[v] += f_scalar[n+v] * x_scalar[n+v];                          \
    ov3 += f_scalar[n+1] * x_scalar[n+0];                               \
    ov4 += f_scalar[n+2] * x_scalar[n+0];                               \
    ov5 += f_scalar[n+2] * x_scalar[n+1];                               \
  }                                                                     \
  ov0 += ovv[0];                                                        \
  ov1 += ovv[1];                                                        \
  ov2 += ovv[2];                                                        \
  if (vwidth > 4) {                                                     \
    ov0 += ovv[4];                                                      \
    ov1 += ovv[5];                                                      \
    ov2 += ovv[6];                                                      \
  }                                                                     \
  if (vwidth > 8) {                                                     \
    ov0 += ovv[8] + ovv[12];                                            \
    ov1 += ovv[9] + ovv[13];                                            \
    ov2 += ovv[10] + ovv[14];                                           \
  }                                                                     \
}

#define IP_PRE_fdotr_acc_force(nall, minlocal, nthreads, f_start,       \
                               f_stride, pos, offload, vflag, ov0, ov1, \
                               ov2, ov3, ov4, ov5)                      \
{                                                                       \
  int o_range = (nall - minlocal) * 4;                                  \
  IP_PRE_omp_range_id_align(iifrom, iito, tid, o_range, nthreads,       \
                            sizeof(acc_t));                             \
                                                                        \
  acc_t *f_scalar = &f_start[0].x;                                      \
  int f_stride4 = f_stride * 4;                                         \
  int t;                                                                \
  if (vflag == VIRIAL_FDOTR) t = 4; else t = 1;                         \
  acc_t *f_scalar2 = f_scalar + f_stride4 * t;                          \
  for ( ; t < nthreads; t++) {                                          \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("simd")                                            \
    for (int n = iifrom; n < iito; n++)                                 \
      f_scalar[n] += f_scalar2[n];                                      \
    f_scalar2 += f_stride4;                                             \
  }                                                                     \
                                                                        \
  if (vflag == VIRIAL_FDOTR) {                                          \
    int nt_min = MIN(4,nthreads);                                       \
    IP_PRE_fdotr_acc_force_l5(iifrom, iito, minlocal, nt_min, f_start,  \
                              f_stride, pos, ov0, ov1, ov2, ov3, ov4,   \
                              ov5);                                     \
  }                                                                     \
}

#else

#define IP_PRE_fdotr_acc_force_l5(lf, lt, minlocal, nthreads, f_start,  \
                                  f_stride, pos, ov0, ov1, ov2,         \
                                  ov3, ov4, ov5)                        \
{                                                                       \
  acc_t *f_scalar = &f_start[0].x;                                      \
  flt_t *x_scalar = &pos[minlocal].x;                                   \
  int f_stride4 = f_stride * 4;                                         \
  _alignvar(acc_t ovv[16],64);                                          \
  int vwidth;                                                           \
  if (sizeof(acc_t) == sizeof(double))                                  \
    vwidth = INTEL_COMPILE_WIDTH/2;                                     \
  else                                                                  \
    vwidth = INTEL_COMPILE_WIDTH;                                       \
  if (vwidth < 4) vwidth = 4;                                           \
  _use_simd_pragma("omp simd aligned(ovv:64)")                          \
  for (int v = 0; v < vwidth; v++) ovv[v] = (acc_t)0.0;                 \
  int remainder = lt % vwidth;                                          \
  if (lf > lt) remainder = 0;                                           \
  const int v_range = lt - remainder;                                   \
  if (nthreads == 2) {                                                  \
    acc_t *f_scalar2 = f_scalar + f_stride4;                            \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("omp simd aligned(f_scalar,f_scalar2,ovv,x_scalar:64)")\
      for (int v = 0; v < vwidth; v++) {                                \
        f_scalar[n+v] += f_scalar2[n+v];                                \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      }                                                                 \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    _use_simd_pragma("loop_count min(4) max(INTEL_COMPILE_WIDTH)")      \
    for (int n = v_range; n < lt; n++)                                  \
      f_scalar[n] += f_scalar2[n];                                      \
  } else if (nthreads==4) {                                             \
    acc_t *f_scalar2 = f_scalar + f_stride4;                            \
    acc_t *f_scalar3 = f_scalar2 + f_stride4;                           \
    acc_t *f_scalar4 = f_scalar3 + f_stride4;                           \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("omp simd aligned(f_scalar,f_scalar2,f_scalar3,f_scalar4,ovv:64)") \
      for (int v = 0; v < vwidth; v++) {                                \
        f_scalar[n+v] += f_scalar2[n+v] + f_scalar3[n+v] +              \
          f_scalar4[n+v];                                               \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      }                                                                 \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    _use_simd_pragma("loop_count min(4) max(INTEL_COMPILE_WIDTH)")      \
    for (int n = v_range; n < lt; n++)                                  \
      f_scalar[n] += f_scalar2[n] + f_scalar3[n] + f_scalar4[n];        \
  } else if (nthreads==1) {                                             \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("omp simd aligned(ovv,f_scalar,x_scalar:64)")    \
      for (int v = 0; v < vwidth; v++)                                  \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
  } else if (nthreads==3) {                                             \
    acc_t *f_scalar2 = f_scalar + f_stride4;                            \
    acc_t *f_scalar3 = f_scalar2 + f_stride4;                           \
    for (int n = lf; n < v_range; n += vwidth) {                        \
      _use_simd_pragma("omp simd aligned(f_scalar,f_scalar2,f_scalar3,ovv,x_scalar:64)") \
      for (int v = 0; v < vwidth; v++) {                                \
        f_scalar[n+v] += f_scalar2[n+v] + f_scalar3[n+v];               \
        ovv[v] += f_scalar[n+v] * x_scalar[n+v];                        \
      }                                                                 \
      ov3 += f_scalar[n+1] * x_scalar[n+0];                             \
      ov4 += f_scalar[n+2] * x_scalar[n+0];                             \
      ov5 += f_scalar[n+2] * x_scalar[n+1];                             \
      if (vwidth > 4) {                                                 \
        ov3 += f_scalar[n+5] * x_scalar[n+4];                           \
        ov4 += f_scalar[n+6] * x_scalar[n+4];                           \
        ov5 += f_scalar[n+6] * x_scalar[n+5];                           \
      }                                                                 \
      if (vwidth > 8) {                                                 \
        ov3 += f_scalar[n+9] * x_scalar[n+8];                           \
        ov3 += f_scalar[n+13] * x_scalar[n+12];                         \
        ov4 += f_scalar[n+10] * x_scalar[n+8];                          \
        ov4 += f_scalar[n+14] * x_scalar[n+12];                         \
        ov5 += f_scalar[n+10] * x_scalar[n+9];                          \
        ov5 += f_scalar[n+14] * x_scalar[n+13];                         \
      }                                                                 \
    }                                                                   \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    _use_simd_pragma("loop_count min(4) max(INTEL_COMPILE_WIDTH)")      \
    for (int n = v_range; n < lt; n++)                                  \
      f_scalar[n] += f_scalar2[n] + f_scalar3[n];                       \
  }                                                                     \
  for (int n = v_range; n < lt; n += 4) {                               \
    _use_simd_pragma("vector aligned")                                  \
    _use_simd_pragma("ivdep")                                           \
    for (int v = 0; v < 4; v++)                                         \
      ovv[v] += f_scalar[n+v] * x_scalar[n+v];                          \
    ov3 += f_scalar[n+1] * x_scalar[n+0];                               \
    ov4 += f_scalar[n+2] * x_scalar[n+0];                               \
    ov5 += f_scalar[n+2] * x_scalar[n+1];                               \
  }                                                                     \
  ov0 += ovv[0];                                                        \
  ov1 += ovv[1];                                                        \
  ov2 += ovv[2];                                                        \
  if (vwidth > 4) {                                                     \
    ov0 += ovv[4];                                                      \
    ov1 += ovv[5];                                                      \
    ov2 += ovv[6];                                                      \
  }                                                                     \
  if (vwidth > 8) {                                                     \
    ov0 += ovv[8] + ovv[12];                                            \
    ov1 += ovv[9] + ovv[13];                                            \
    ov2 += ovv[10] + ovv[14];                                           \
  }                                                                     \
}

#define IP_PRE_fdotr_acc_force(nall, minlocal, nthreads, f_start,       \
                               f_stride, pos, offload, vflag, ov0, ov1, \
                               ov2, ov3, ov4, ov5)                      \
{                                                                       \
  int o_range = (nall - minlocal) * 4;                                  \
  IP_PRE_omp_range_id_align(iifrom, iito, tid, o_range, nthreads,       \
                            sizeof(acc_t));                             \
                                                                        \
  acc_t *f_scalar = &f_start[0].x;                                      \
  int f_stride4 = f_stride * 4;                                         \
  int t;                                                                \
  if (vflag == VIRIAL_FDOTR) t = 4; else t = 1;                         \
  acc_t *f_scalar2 = f_scalar + f_stride4 * t;                          \
  for ( ; t < nthreads; t++) {                                          \
    _use_simd_pragma("omp simd aligned(f_scalar,f_scalar2:64)")         \
    for (int n = iifrom; n < iito; n++)                                 \
      f_scalar[n] += f_scalar2[n];                                      \
    f_scalar2 += f_stride4;                                             \
  }                                                                     \
                                                                        \
  if (vflag == VIRIAL_FDOTR) {                                          \
    int nt_min = MIN(4,nthreads);                                       \
    IP_PRE_fdotr_acc_force_l5(iifrom, iito, minlocal, nt_min, f_start,  \
                              f_stride, pos, ov0, ov1, ov2, ov3, ov4,   \
                              ov5);                                     \
  }                                                                     \
}

#endif

#ifdef _LMP_INTEL_OFFLOAD
#include <sys/time.h>

__declspec( target (mic))
inline double MIC_Wtime() {
  double time;
  struct timeval tv;

  gettimeofday(&tv, nullptr);
  time = 1.0 * tv.tv_sec + 1.0e-6 * tv.tv_usec;
  return time;
}

#define IP_PRE_neighbor_pad(jnum, offload)                              \
{                                                                       \
  const int opad_mask = ~static_cast<int>(INTEL_MIC_NBOR_PAD *          \
                                          sizeof(float) /               \
                                          sizeof(flt_t) - 1);           \
  const int pad_mask = ~static_cast<int>(INTEL_NBOR_PAD *               \
                                          sizeof(float) /               \
                                          sizeof(flt_t) - 1);           \
  if (offload && INTEL_MIC_NBOR_PAD > 1)                                \
    jnum = (jnum + INTEL_MIC_NBOR_PAD * sizeof(float) /                 \
            sizeof(flt_t) - 1) & opad_mask;                             \
  else if (INTEL_NBOR_PAD > 1)                                          \
    jnum = (jnum + INTEL_NBOR_PAD * sizeof(float) /                     \
            sizeof(flt_t) - 1) & pad_mask;                              \
}

#define IP_PRE_pack_separate_buffers(fix, buffers, ago, offload,        \
                                     nlocal, nall)                      \
{                                                                       \
    if (fix->separate_buffers() && ago != 0) {                          \
    fix->start_watch(TIME_PACK);                                        \
    if (offload) {                                                      \
      int packthreads;                                                  \
      if (comm->nthreads > INTEL_HTHREADS) packthreads = comm->nthreads;\
      else packthreads = 1;                                             \
      _use_omp_pragma("omp parallel if(packthreads > 1)")               \
      {                                                                 \
        int ifrom, ito, tid;                                            \
        IP_PRE_omp_range_id_align(ifrom, ito, tid, nlocal,              \
                                  packthreads, sizeof(flt_t));          \
        buffers->thr_pack_cop(ifrom, ito, 0);                           \
        int nghost = nall - nlocal;                                     \
        if (nghost) {                                                   \
          IP_PRE_omp_range_align(ifrom, ito, tid, nall - nlocal,        \
                                 packthreads, sizeof(flt_t));           \
          buffers->thr_pack_cop(ifrom + nlocal, ito + nlocal,           \
                                fix->offload_min_ghost() - nlocal,      \
                                ago == 1);                              \
        }                                                               \
      }                                                                 \
    } else {                                                            \
      buffers->thr_pack_host(fix->host_min_local(), nlocal, 0);         \
      buffers->thr_pack_host(nlocal, nall,                              \
                             fix->host_min_ghost()-nlocal);             \
    }                                                                   \
    fix->stop_watch(TIME_PACK);                                         \
  }                                                                     \
}

#define IP_PRE_get_transfern(ago, newton, eflag, vflag,                 \
                             buffers, offload, fix, separate_flag,      \
                             x_size, q_size, ev_size, f_stride)         \
{                                                                       \
  separate_flag = 0;                                                    \
  if (ago == 0) {                                                       \
    x_size = 0;                                                         \
    q_size = nall;                                                      \
    if (offload) {                                                      \
      if (fix->separate_buffers()) {                                    \
        if (lmp->atom->torque)                                          \
          separate_flag = 2;                                            \
        else                                                            \
          separate_flag = 1;                                            \
      } else                                                            \
        separate_flag = 3;                                              \
    }                                                                   \
  } else {                                                              \
    x_size = nall;                                                      \
    q_size = 0;                                                         \
  }                                                                     \
  ev_size = 0;                                                          \
  if (eflag) ev_size = 2;                                               \
  if (vflag) ev_size = 8;                                               \
  if (newton)                                                           \
    f_stride = buffers->get_stride(nall);                               \
  else                                                                  \
    f_stride = buffers->get_stride(nlocal);                             \
}

#define IP_PRE_get_buffers(offload, buffers, fix, tc, f_start,          \
                           ev_global)                                   \
{                                                                       \
  if (offload) {                                                        \
    tc = buffers->get_off_threads();                                    \
    f_start = buffers->get_off_f();                                     \
    ev_global = buffers->get_ev_global();                               \
  } else {                                                              \
    tc = comm->nthreads;                                                \
    f_start = buffers->get_f();                                         \
    fix->start_watch(TIME_HOST_PAIR);                                   \
    ev_global = buffers->get_ev_global_host();                          \
  }                                                                     \
}

#define IP_PRE_repack_for_offload(newton, separate_flag, nlocal, nall,  \
                                  f_stride, x, q)                       \
{                                                                       \
  if (separate_flag) {                                                  \
    if (separate_flag < 3) {                                            \
      int all_local = nlocal;                                           \
      int ghost_min = overflow[LMP_GHOST_MIN];                          \
      nlocal = overflow[LMP_LOCAL_MAX] + 1;                             \
      int nghost = overflow[LMP_GHOST_MAX] + 1 - ghost_min;             \
      if (nghost < 0) nghost = 0;                                       \
      nall = nlocal + nghost;                                           \
      separate_flag--;                                                  \
      int flength;                                                      \
      if (newton) flength = nall;                                       \
      else flength = nlocal;                                            \
      IP_PRE_get_stride(f_stride, flength, sizeof(FORCE_T),             \
                           separate_flag);                              \
      if (nghost) {                                                     \
        if (nlocal < all_local || ghost_min > all_local) {              \
          memmove(x + nlocal, x + ghost_min,                            \
                  (nall - nlocal) * sizeof(ATOM_T));                    \
          if (q != 0)                                                   \
            memmove((void *)(q + nlocal), (void *)(q + ghost_min),      \
                    (nall - nlocal) * sizeof(flt_t));                   \
        }                                                               \
      }                                                                 \
    }                                                                   \
    x[nall].x = INTEL_BIGP;                                             \
    x[nall].y = INTEL_BIGP;                                             \
    x[nall].z = INTEL_BIGP;                                             \
  }                                                                     \
}

#define IP_PRE_fdotr_reduce_omp(newton, nall, minlocal, nthreads,       \
                                f_start, f_stride, x, offload, vflag,   \
                                ov0, ov1, ov2, ov3, ov4, ov5)           \
{                                                                       \
  if (newton) {                                                         \
    _use_omp_pragma("omp barrier");                                     \
    IP_PRE_fdotr_acc_force(nall, minlocal, nthreads, f_start,           \
                           f_stride, x, offload, vflag, ov0, ov1, ov2,  \
                           ov3, ov4, ov5);                              \
  }                                                                     \
}

#define IP_PRE_fdotr_reduce(newton, nall, nthreads, f_stride, vflag,    \
                            ov0, ov1, ov2, ov3, ov4, ov5)

#else

#if INTEL_NBOR_PAD > 1

#define IP_PRE_neighbor_pad(jnum, offload)                              \
{                                                                       \
  const int pad_mask = ~static_cast<int>(INTEL_NBOR_PAD *               \
                                         sizeof(float) /                \
                                         sizeof(flt_t) - 1);            \
  jnum = (jnum + INTEL_NBOR_PAD * sizeof(float) /                       \
          sizeof(flt_t) - 1) & pad_mask;                                \
}

#else

#define IP_PRE_neighbor_pad(jnum, offload)

#endif

#define MIC_Wtime MPI_Wtime
#define IP_PRE_pack_separate_buffers(fix, buffers, ago, offload,        \
                                     nlocal, nall)

#define IP_PRE_get_transfern(ago, newton, eflag, vflag,                 \
                             buffers, offload, fix, separate_flag,      \
                             x_size, q_size, ev_size, f_stride)         \
{                                                                       \
  separate_flag = 0;                                                    \
  int f_length;                                                         \
  if (newton)                                                           \
    f_length = nall;                                                    \
  else                                                                  \
    f_length = nlocal;                                                  \
  f_stride = buffers->get_stride(f_length);                             \
}

#define IP_PRE_get_buffers(offload, buffers, fix, tc, f_start,          \
                           ev_global)                                   \
{                                                                       \
  tc = comm->nthreads;                                                  \
  f_start = buffers->get_f();                                           \
  fix->start_watch(TIME_HOST_PAIR);                                     \
  ev_global = buffers->get_ev_global_host();                            \
}

#define IP_PRE_repack_for_offload(newton, separate_flag, nlocal, nall,  \
                                  f_stride, x, q)

#define IP_PRE_fdotr_reduce_omp(newton, nall, minlocal, nthreads,       \
                                f_start, f_stride, x, offload, vflag,   \
                                ov0, ov1, ov2, ov3, ov4, ov5)           \
{                                                                       \
  if (newton) {                                                         \
    if (vflag == 2 && nthreads > INTEL_HTHREADS) {                      \
      _use_omp_pragma("omp barrier");                                   \
      buffers->fdotr_reduce(nall, nthreads, f_stride, ov0, ov1, ov2,    \
                            ov3, ov4, ov5);                             \
    }                                                                   \
  }                                                                     \
}

#define IP_PRE_fdotr_reduce(newton, nall, nthreads, f_stride, vflag,    \
                            ov0, ov1, ov2, ov3, ov4, ov5)               \
{                                                                       \
  if (newton) {                                                         \
    if (vflag == 2 && nthreads <= INTEL_HTHREADS) {                     \
      int lt = nall * 4;                                                \
      buffers->fdotr_reduce_l5(0, lt, nthreads, f_stride, ov0, ov1,     \
                               ov2, ov3, ov4, ov5);                     \
    }                                                                   \
  }                                                                     \
}

#endif

#define IP_PRE_ev_tally_nbor(vflag, fpair, delx, dely, delz)            \
{                                                                       \
  if (vflag == 1) {                                                     \
    sv0 += delx * delx * fpair;                                         \
    sv1 += dely * dely * fpair;                                         \
    sv2 += delz * delz * fpair;                                         \
    sv3 += delx * dely * fpair;                                         \
    sv4 += delx * delz * fpair;                                         \
    sv5 += dely * delz * fpair;                                         \
  }                                                                     \
}

#define IP_PRE_ev_tally_nborv(vflag, dx, dy, dz, fpx, fpy, fpz)         \
{                                                                       \
  if (vflag == 1) {                                                     \
    sv0 += dx * fpx;                                                    \
    sv1 += dy * fpy;                                                    \
    sv2 += dz * fpz;                                                    \
    sv3 += dx * fpy;                                                    \
    sv4 += dx * fpz;                                                    \
    sv5 += dy * fpz;                                                    \
  }                                                                     \
}

#define IP_PRE_ev_tally_nbor3(vflag, fj, fk, delx, dely, delz, delr2)   \
{                                                                       \
  if (vflag == 1) {                                                     \
    sv0 += delx * fj[0] + delr2[0] * fk[0];                             \
    sv1 += dely * fj[1] + delr2[1] * fk[1];                             \
    sv2 += delz * fj[2] + delr2[2] * fk[2];                             \
    sv3 += delx * fj[1] + delr2[0] * fk[1];                             \
    sv4 += delx * fj[2] + delr2[0] * fk[2];                             \
    sv5 += dely * fj[2] + delr2[1] * fk[2];                             \
  }                                                                     \
}

#define IP_PRE_ev_tally_nbor3v(vflag, fj0, fj1, fj2, delx, dely, delz)  \
{                                                                       \
  if (vflag == 1) {                                                     \
    sv0 += delx * fj0;                                                  \
    sv1 += dely * fj1;                                                  \
    sv2 += delz * fj2;                                                  \
    sv3 += delx * fj1;                                                  \
    sv4 += delx * fj2;                                                  \
    sv5 += dely * fj2;                                                  \
  }                                                                     \
}

#define IP_PRE_ev_tally_bond(eflag, VFLAG, eatom, vflag, ebond, i1, i2, \
                             fbond, delx, dely, delz, obond, force,     \
                             newton, nlocal, ov0, ov1, ov2, ov3, ov4,   \
                             ov5)                                       \
{                                                                       \
  flt_t ev_pre;                                                         \
  if (newton) ev_pre = (flt_t)1.0;                                      \
  else {                                                                \
    ev_pre = (flt_t)0.0;                                                \
    if (i1 < nlocal) ev_pre += (flt_t)0.5;                              \
    if (i2 < nlocal) ev_pre += (flt_t)0.5;                              \
  }                                                                     \
                                                                        \
  if (eflag) {                                                          \
    obond += ev_pre * ebond;                                            \
    if (eatom) {                                                        \
      flt_t halfeng = ebond * (flt_t)0.5;                               \
      if (newton || i1 < nlocal) f[i1].w += halfeng;                    \
      if (newton || i2 < nlocal) f[i2].w += halfeng;                    \
    }                                                                   \
  }                                                                     \
                                                                        \
  if (VFLAG && vflag) {                                                 \
    ov0 += ev_pre * (delx * delx * fbond);                              \
    ov1 += ev_pre * (dely * dely * fbond);                              \
    ov2 += ev_pre * (delz * delz * fbond);                              \
    ov3 += ev_pre * (delx * dely * fbond);                              \
    ov4 += ev_pre * (delx * delz * fbond);                              \
    ov5 += ev_pre * (dely * delz * fbond);                              \
  }                                                                     \
}

#define IP_PRE_ev_tally_angle(eflag, VFLAG, eatom, vflag, eangle, i1,   \
                              i2, i3, f1x, f1y, f1z, f3x, f3y, f3z,     \
                              delx1, dely1, delz1, delx2, dely2, delz2, \
                              oeangle, force, newton, nlocal, ov0, ov1, \
                              ov2, ov3, ov4, ov5)                       \
{                                                                       \
  flt_t ev_pre;                                                         \
  if (newton) ev_pre = (flt_t)1.0;                                      \
  else {                                                                \
    ev_pre = (flt_t)0.0;                                                \
    if (i1 < nlocal) ev_pre += (flt_t)0.3333333333333333;               \
    if (i2 < nlocal) ev_pre += (flt_t)0.3333333333333333;               \
    if (i3 < nlocal) ev_pre += (flt_t)0.3333333333333333;               \
  }                                                                     \
                                                                        \
  if (eflag) {                                                          \
    oeangle += ev_pre * eangle;                                         \
    if (eatom) {                                                        \
      flt_t thirdeng = eangle * (flt_t)0.3333333333333333;              \
      if (newton || i1 < nlocal) f[i1].w += thirdeng;                   \
      if (newton || i2 < nlocal) f[i2].w += thirdeng;                   \
      if (newton || i3 < nlocal) f[i3].w += thirdeng;                   \
    }                                                                   \
  }                                                                     \
                                                                        \
  if (VFLAG && vflag) {                                                 \
    ov0 += ev_pre * (delx1 * f1x + delx2 * f3x);                        \
    ov1 += ev_pre * (dely1 * f1y + dely2 * f3y);                        \
    ov2 += ev_pre * (delz1 * f1z + delz2 * f3z);                        \
    ov3 += ev_pre * (delx1 * f1y + delx2 * f3y);                        \
    ov4 += ev_pre * (delx1 * f1z + delx2 * f3z);                        \
    ov5 += ev_pre * (dely1 * f1z + dely2 * f3z);                        \
  }                                                                     \
}

#define IP_PRE_ev_tally_dihed(eflag, VFLAG, eatom, vflag, deng, i1, i2, \
                              i3, i4, f1x, f1y, f1z, f3x, f3y, f3z, f4x,\
                              f4y, f4z, vb1x, vb1y, vb1z, vb2x, vb2y,   \
                              vb2z, vb3x, vb3y, vb3z, oedihedral, force,\
                              newton, nlocal, ov0, ov1, ov2, ov3, ov4,  \
                              ov5)                                      \
{                                                                       \
  flt_t ev_pre;                                                         \
  if (newton) ev_pre = (flt_t)1.0;                                      \
  else {                                                                \
    ev_pre = (flt_t)0.0;                                                \
    if (i1 < nlocal) ev_pre += (flt_t)0.25;                             \
    if (i2 < nlocal) ev_pre += (flt_t)0.25;                             \
    if (i3 < nlocal) ev_pre += (flt_t)0.25;                             \
    if (i4 < nlocal) ev_pre += (flt_t)0.25;                             \
  }                                                                     \
                                                                        \
  if (eflag) {                                                          \
    oedihedral += ev_pre * deng;                                        \
    if (eatom) {                                                        \
      flt_t qdeng = deng * (flt_t)0.25;                                 \
      if (newton || i1 < nlocal) f[i1].w += qdeng;                      \
      if (newton || i2 < nlocal) f[i2].w += qdeng;                      \
      if (newton || i3 < nlocal) f[i3].w += qdeng;                      \
      if (newton || i4 < nlocal) f[i4].w += qdeng;                      \
    }                                                                   \
  }                                                                     \
                                                                        \
  if (VFLAG && vflag) {                                                 \
    ov0 += ev_pre * (vb1x*f1x + vb2x*f3x + (vb3x+vb2x)*f4x);            \
    ov1 += ev_pre * (vb1y*f1y + vb2y*f3y + (vb3y+vb2y)*f4y);            \
    ov2 += ev_pre * (vb1z*f1z + vb2z*f3z + (vb3z+vb2z)*f4z);            \
    ov3 += ev_pre * (vb1x*f1y + vb2x*f3y + (vb3x+vb2x)*f4y);            \
    ov4 += ev_pre * (vb1x*f1z + vb2x*f3z + (vb3x+vb2x)*f4z);            \
    ov5 += ev_pre * (vb1y*f1z + vb2y*f3z + (vb3y+vb2y)*f4z);            \
  }                                                                     \
}

#define IP_PRE_ev_tally_atom(newton, eflag, vflag, f, fwtmp)            \
{                                                                       \
  if (eflag) {                                                          \
    f[i].w += fwtmp;                                                    \
    oevdwl += sevdwl;                                                   \
  }                                                                     \
  if (newton == 0 && vflag == 1) {                                      \
    ov0 += sv0;                                                         \
    ov1 += sv1;                                                         \
    ov2 += sv2;                                                         \
    ov3 += sv3;                                                         \
    ov4 += sv4;                                                         \
    ov5 += sv5;                                                         \
  }                                                                     \
}

#define IP_PRE_ev_tally_atomq(newton, eflag, vflag, f, fwtmp)           \
{                                                                       \
  if (eflag) {                                                          \
    f[i].w += fwtmp;                                                    \
    oevdwl += sevdwl;                                                   \
    oecoul += secoul;                                                   \
  }                                                                     \
  if (newton == 0 && vflag == 1) {                                      \
    ov0 += sv0;                                                         \
    ov1 += sv1;                                                         \
    ov2 += sv2;                                                         \
    ov3 += sv3;                                                         \
    ov4 += sv4;                                                         \
    ov5 += sv5;                                                         \
  }                                                                     \
}

}

#endif
