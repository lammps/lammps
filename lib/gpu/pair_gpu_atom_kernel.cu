/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
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

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp4 double4
#else
#define numtyp float
#define numtyp4 float4
#endif

#ifdef NV_KERNEL
#include "geryon/ucl_nv_kernel.h"
#else
#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define GLOBAL_ID_X get_global_id(0)
#endif

__kernel void kernel_cast_x(__global numtyp4 *x_type, __global double *x,
                            __global int *type, const int nall) {
  int ii=GLOBAL_ID_X;

  if (ii<nall) {
    numtyp4 xt;
    xt.w=type[ii];
    int i=ii*3;
    xt.x=x[i];
    xt.y=x[i+1];
    xt.z=x[i+2];
    x_type[ii]=xt;
  } // if ii
}
