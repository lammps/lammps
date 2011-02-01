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

#ifndef PPPM_GPU_KERNEL
#define PPPM_GPU_KERNEL

#define OFFSET 16384

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp2 double2
#define numtyp4 double4
#define acctyp double
#define acctyp4 double4
#endif

#ifdef _SINGLE_DOUBLE
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp double
#define acctyp4 double4
#endif

#ifndef numtyp
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp float
#define acctyp4 float4
#endif

#ifdef NV_KERNEL

#include "geryon/ucl_nv_kernel.h"
texture<float4> pos_tex;
texture<float> q_tex;

#ifdef _DOUBLE_DOUBLE
__inline double4 fetch_pos(const int& i, const double4 *pos)
{
  return pos[i];
}
__inline double fetch_q(const int& i, const double *q)
{
  return q[i];
}
#else
__inline float4 fetch_pos(const int& i, const float4 *pos)
{
  return tex1Dfetch(pos_tex, i);
}
__inline float fetch_q(const int& i, const float *q)
{
  return tex1Dfetch(q_tex, i);
}
#endif

#else

#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define GLOBAL_ID_X get_global_id(0)
#define THREAD_ID_X get_local_id(0)
#define BLOCK_ID_X get_group_id(0)
#define BLOCK_SIZE_X get_local_size(0)
#define __syncthreads() barrier(CLK_LOCAL_MEM_FENCE)
#define __inline inline

#define fetch_pos(i,y) x_[i]
#define fetch_q(i,y) q_[i]

#endif

__kernel void particle_map(__global numtyp4 *x_, const int nlocal, 
                           __global int *counts, __global int *ans, 
                           const numtyp boxlo_x, const numtyp boxlo_y,
                           const numtyp boxlo_z, const numtyp delxinv,
                           const numtyp delyinv, const numtyp delzinv,
                           const numtyp shift, const int nxlo_out,
                           const int nxhi_out, const int nylo_out,
                           const int nyhi_out, const int nzlo_out,
                           const int nzhi_out, const int nlower,
                           const int nupper, __global int *error) {
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;
  int nx,ny,nz;
/*
  if (ii<nlocal) {
    numtyp4 p=fetch_pos(ii,x_);
// shift boxlo

    nx = int((p.x-boxlo_x)*delxinv+shift) - OFFSET;
    ny = int((p.y-boxlo_y)*delyinv+shift) - OFFSET;
    nz = int((p.z-boxlo_z)*delzinv+shift) - OFFSET;
    counts[
    part2grid[i][0] = nx;
    part2grid[i][1] = ny;
    part2grid[i][2] = nz;

    // check that entire stencil around nx,ny,nz will fit in my 3d brick
    if (nx+nlower < nxlo_out || nx+nupper > nxhi_out ||
        ny+nlower < nylo_out || ny+nupper > nyhi_out ||
        nz+nlower < nzlo_out || nz+nupper > nzhi_out) *error=1;
  }*/
}

#endif

