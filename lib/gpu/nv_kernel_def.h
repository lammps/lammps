/* ----------------------------------------------------------------------
   LAMMPS-Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

/*************************************************************************
                 See pair_gpu_dev_kernel.cu for definitions
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
