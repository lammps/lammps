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
                           Preprocessor Definitions
                           
  Note: It is assumed that constants with the same names are defined with
  the same values in all files.
  
  ARCH
     Definition:   Architecture number for accelerator
  MEM_THREADS
     Definition:   Number of threads with sequential ids accessing memory
                   simultaneously on multiprocessor
  WARP_SIZE:
     Definition:   Number of threads guaranteed to be on the same instruction
  THREADS_PER_ATOM
     Definition:   Default number of threads assigned per atom for pair styles
     Restructions: Must be power of 2; THREADS_PER_ATOM<=WARP_SIZE
  THREADS_PER_CHARGE
     Definition:   Default number of threads assigned per atom for pair styles
                   with charge
     Restructions: Must be power of 2; THREADS_PER_ATOM<=WARP_SIZE
  PPPM_MAX_SPLINE
     Definition:   Maximum order for splines in PPPM
  PPPM_BLOCK_1D    
     Definition:   Thread block size for PPPM kernels
     Restrictions: PPPM_BLOCK_1D>=PPPM_MAX_SPLINE*PPPM_MAX_SPLINE
                   PPPM_BLOCK_1D%32==0 
  BLOCK_PAIR
     Definition:   Default thread block size for pair styles
     Restrictions:
  MAX_SHARED_TYPES 8
     Definition:   Max number of atom type params can be stored in shared memory
     Restrictions: MAX_SHARED_TYPES*MAX_SHARED_TYPES<=BLOCK_PAIR
  BLOCK_CELL_2D 
     Definition:   Default block size in each dimension for cell list builds
                   and matrix transpose
  BLOCK_CELL_ID    
     Definition:   Default block size for binning atoms in cell list builds
  BLOCK_NBOR_BUILD 
     Definition:   Default block size for neighbor list builds
  BLOCK_BIO_PAIR
     Definition:   Default thread block size for "bio" pair styles
  MAX_BIO_SHARED_TYPES
     Definition:   Max number of atom type params can be stored in shared memory
     Restrictions:  MAX_BIO_SHARED_TYPES<=BLOCK_BIO_PAIR*2 &&
                    MAX_BIO_SHARED_TYPES>=BLOCK_BIO_PAIR

*************************************************************************/

#ifndef PAIR_GPU_DEV_KERNEL
#define PAIR_GPU_DEV_KERNEL

#ifdef NV_KERNEL

#include "nv_kernel_def.h"

#else

#define GLOBAL_ID_X get_global_id(0)
#define ARCH 0
#define DRIVER 0
#define MEM_THREADS 16
#define WARP_SIZE 1
#define THREADS_PER_ATOM 1
#define THREADS_PER_CHARGE 1
#define BLOCK_PAIR 64
#define MAX_SHARED_TYPES 8
#define BLOCK_NBOR_BUILD 64
#define BLOCK_BIO_PAIR 64

#endif

#define PPPM_MAX_SPLINE 8
#define PPPM_BLOCK_1D 64
#define BLOCK_CELL_2D 8
#define BLOCK_CELL_ID 128
#define MAX_BIO_SHARED_TYPES 128

__kernel void kernel_zero(__global int *mem, int numel) {
  int ii=GLOBAL_ID_X;
  
  if (ii<numel)
    mem[ii]=0;
}

__kernel void kernel_info(__global int *info) {
  info[0]=ARCH;
  info[1]=MEM_THREADS;
  info[2]=WARP_SIZE;
  info[3]=THREADS_PER_ATOM;
  info[4]=PPPM_MAX_SPLINE;
  info[5]=PPPM_BLOCK_1D;
  info[6]=BLOCK_PAIR;
  info[7]=MAX_SHARED_TYPES;
  info[8]=BLOCK_CELL_2D;
  info[9]=BLOCK_CELL_ID;
  info[10]=BLOCK_NBOR_BUILD;
  info[11]=BLOCK_BIO_PAIR;
  info[12]=MAX_BIO_SHARED_TYPES;
  info[13]=THREADS_PER_CHARGE;
}

#endif

