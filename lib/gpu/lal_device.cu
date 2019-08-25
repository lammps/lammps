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

#ifdef NV_KERNEL
#include "lal_preprocessor.h"
#endif

__kernel void kernel_zero(__global int *restrict mem,
                          int numel) {
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
  info[14]=BLOCK_ELLIPSE;
}

