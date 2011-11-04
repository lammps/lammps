/***************************************************************************
                               ucl_nv_kernel.h
                             -------------------
                               W. Michael Brown

  Preprocessor macros for OpenCL/CUDA compatibility

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Mon May 3 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

// Only allow this file to be included by CUDA and OpenCL specific headers
#ifndef UCL_NV_KERNEL_H
#define UCL_NV_KERNEL_H

#if (__CUDA_ARCH__ < 200)
#define mul24 __mul24
#define MEM_THREADS 16
#else
#define mul24(X,Y) (X)*(Y)
#define MEM_THREADS 32
#endif

#ifdef CUDA_PRE_THREE
struct __builtin_align__(16) _double4
{
  double x, y, z, w;
};
typedef struct _double4 double4;
#endif

#define GLOBAL_ID_X threadIdx.x+mul24(blockIdx.x,blockDim.x)
#define GLOBAL_ID_Y threadIdx.y+mul24(blockIdx.y,blockDim.y)
#define GLOBAL_SIZE_X mul24(gridDim.x,blockDim.x);
#define GLOBAL_SIZE_Y mul24(gridDim.y,blockDim.y);
#define THREAD_ID_X threadIdx.x
#define THREAD_ID_Y threadIdx.y
#define BLOCK_ID_X blockIdx.x
#define BLOCK_ID_Y blockIdx.y
#define BLOCK_SIZE_X blockDim.x
#define BLOCK_SIZE_Y blockDim.y
#define __kernel extern "C" __global__
#define __local __shared__
#define __global  
#define atom_add atomicAdd
#define ucl_inline static __inline__ __device__ 

#endif

