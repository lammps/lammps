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

#ifndef PAIR_GPU_DEV_KERNEL
#define PAIR_GPU_DEV_KERNEL

#ifdef NV_KERNEL

#include "geryon/ucl_nv_kernel.h"
#ifdef __CUDA_ARCH__
#define ARCH __CUDA_ARCH__
#else
#define ARCH 100
#endif

#else

#define GLOBAL_ID_X get_global_id(0)
#define ARCH 0
#define DRIVER 0
#define MEM_THREADS 16;

#endif

__kernel void kernel_zero(__global int *mem, int numel) {
  int ii=GLOBAL_ID_X;
  
  if (ii<numel)
    mem[ii]=0;
}

__kernel void kernel_info(__global int *info) {
  info[0]=ARCH;
  info[1]=MEM_THREADS;
}

#endif

