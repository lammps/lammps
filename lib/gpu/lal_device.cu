// **************************************************************************
//                                  device.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for device information
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_preprocessor.h"
#endif

__kernel void kernel_zero(__global int *restrict mem,
                          int numel) {
  int ii=GLOBAL_ID_X;

  if (ii<numel)
    mem[ii]=0;
}

__kernel void kernel_info(__global int *info) {
  #ifdef __CUDA_ARCH__
  info[0]=__CUDA_ARCH__;
  #else
  info[0]=0;
  #endif
  info[1]=CONFIG_ID;
  info[2]=SIMD_SIZE;
  info[3]=MEM_THREADS;
  info[4]=SHUFFLE_AVAIL;
  info[5]=FAST_MATH;

  info[6]=THREADS_PER_ATOM;
  info[7]=THREADS_PER_CHARGE;
  info[8]=THREADS_PER_THREE;

  info[9]=BLOCK_PAIR;
  info[10]=BLOCK_BIO_PAIR;
  info[11]=BLOCK_ELLIPSE;
  info[12]=PPPM_BLOCK_1D;
  info[13]=BLOCK_NBOR_BUILD;
  info[14]=BLOCK_CELL_2D;
  info[15]=BLOCK_CELL_ID;

  info[16]=MAX_SHARED_TYPES;
  info[17]=MAX_BIO_SHARED_TYPES;
  info[18]=PPPM_MAX_SPLINE;
}
