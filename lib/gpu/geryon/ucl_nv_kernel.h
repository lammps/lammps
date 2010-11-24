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

#define GLOBAL_ID_X threadIdx.x+__mul24(blockIdx.x,blockDim.x)
#define GLOBAL_ID_Y threadIdx.y+__mul24(blockIdx.y,blockDim.y)
#define THREAD_ID_X threadIdx.x
#define THREAD_ID_Y threadIdx.y
#define BLOCK_ID_X blockIdx.x
#define BLOCK_ID_Y blockIdx.y
#define BLOCK_SIZE_X blockDim.x
#define BLOCK_SIZE_Y blockDim.y
#define __kernel extern "C" __global__
#define __local __shared__
#define mul24 __mul24
#define __global  
#define __inline static __inline__ __device__ 

#endif
