// **************************************************************************
//                              nv_kernel_def.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for CUDA-specific preprocessor definitions
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : brownw@ornl.gov
// ***************************************************************************/

/*************************************************************************
                 See device.cu for definitions
                       of preprocessor constants
*************************************************************************/

#ifndef NV_KERNEL_DEF
#define NV_KERNEL_DEF

#include "geryon/ucl_nv_kernel.h"
#ifdef __CUDA_ARCH__
#define ARCH __CUDA_ARCH__
#else
#define ARCH 100
#endif

#if (ARCH < 200)

#define THREADS_PER_ATOM 1
#define THREADS_PER_CHARGE 8
#define BLOCK_NBOR_BUILD 64
#define BLOCK_PAIR 64
#define BLOCK_BIO_PAIR 64
#define MAX_SHARED_TYPES 8

#else

#define THREADS_PER_ATOM 1
#define THREADS_PER_CHARGE 8
#define BLOCK_NBOR_BUILD 128
#define BLOCK_PAIR 128
#define BLOCK_BIO_PAIR 128
#define MAX_SHARED_TYPES 11

#endif

#define WARP_SIZE 32

#endif
