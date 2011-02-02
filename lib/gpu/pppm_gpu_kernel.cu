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
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
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
                           const int npts_x, const int npts_y,
                           const int npts_z, const int _brick_stride,
                           const int max_atoms, __global int *error) {
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;
  int nx,ny,nz;
  numtyp tx,ty,tz;

  if (ii<nlocal) {
    numtyp4 p=fetch_pos(ii,x_);

    // Boxlo is adjusted to include ghost cells so that starting index is 0
    tx=(p.x-boxlo_x)*delxinv;
    nx=int(tx);
    ty=(p.y-boxlo_y)*delyinv;
    ny=int(ty);
    tz=(p.z-boxlo_z)*delzinv;
    nz=int(tz);

    if (tx<0 || ty<0 || tz<0 || nx>=npts_x || ny>=npts_y || nz>=npts_z)
      *error=1;
    else {
      int i=nz*npts_y*npts_x+ny*npts_x+nx;
      int old=atom_add(counts+i, 1);
      if (old==max_atoms)
        *error=2;
      else
        ans[_brick_stride*old+i]=ii;
    }
  }
}

/*
__kernel void particle_map(__global numtyp4 *x_, __global numtyp *q_,
                           __global int *counts, __global int *atoms, 
                           const numtyp boxlo_x, const numtyp boxlo_y,
                           const numtyp boxlo_z, const numtyp delxinv,
                           const numtyp delyinv, const numtyp delzinv,
                           const int npts_x, const int npts_y,
                           const int npts_z, const int _brick_stride,
                           const int max_atoms, __global int *error) {
  // ii indexes the two interacting particles in gi
  int xx=THREAD_ID_X;
  int yy=THREAD_ID_Y;
  int bx=BLOCK_ID_X;
  int by=BLOCK_ID_Y;
  int block_size=BLOCK_SIZE_X;
  
  int max_y=BLOCK_ID_Y*block_size+block_size;
  int max_x=BLOCK_ID_X*block_size+block_size;
  
  __local numtyp4 p;
  __local numtyp q,dx,dy,dz;
  __local int brick_i,count,atom_i;
  
  for (int z=-nlower; z<npts_z-order; z++)
    for (int ny=max_y-block_size; ny<max_y; ny++) {
      if (ny>npts_y)
        break;
      brick_i = z*npts_x*npts_y + ny*npts_x + max_x - block_size;
      for (int nx=max_x-block_size; nx<max_x; nx++) {
        if (nx>npts_x)
          break;
        count=counts[brick_i];
        for (int i=0; i<count; i++) {
          int atom_i=atoms[brick_i+i*_brick_stride];
          p=fetch_pos(x_,atom_i);
          q=fetch_q(q_,atom_i);
          
    dx = nx+shiftone - (x[i][0]-boxlo[0])*delxinv;
    dy = ny+shiftone - (x[i][1]-boxlo[1])*delyinv;
    dz = nz+shiftone - (x[i][2]-boxlo[2])*delzinv;
          
          
        
        
        
        brick_i++;
  
  
  int nx=GLOBAL_ID_X+nlower;
  int ny=GLOBAL_ID_Y+nlower;
  int block_size=BLOCK_SIZE_X;
  
    for (

  // (nx,ny,nz) = global coords of grid pt to "lower left" of charge
  // (dx,dy,dz) = distance to "lower left" grid pt
  // (mx,my,mz) = global coords of moving stencil pt

  double *q = atom->q;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {

    nx = part2grid[i][0];
    ny = part2grid[i][1];
    nz = part2grid[i][2];

    compute_rho1d(dx,dy,dz);

    z0 = delvolinv * q[i];
    for (n = nlower; n <= nupper; n++) {
      mz = n+nz;  // z-index of point being updated
      y0 = z0*rho1d[2][n]; 
      for (m = nlower; m <= nupper; m++) {
        my = m+ny;
        x0 = y0*rho1d[1][m];
        for (l = nlower; l <= nupper; l++) {
          mx = l+nx;
          density_brick[mz][my][mx] += x0*rho1d[0][l];
        }
      }
    }
  }
}

void PPPMGPU::compute_rho1d(double dx, double dy, double dz)
{
  int k,l;

  for (k = (1-order)/2; k <= order/2; k++) {
    rho1d[0][k] = 0.0;
    rho1d[1][k] = 0.0;
    rho1d[2][k] = 0.0;
    for (l = order-1; l >= 0; l--) {
      rho1d[0][k] = rho_coeff[l][k] + rho1d[0][k]*dx;
      rho1d[1][k] = rho_coeff[l][k] + rho1d[1][k]*dy;
      rho1d[2][k] = rho_coeff[l][k] + rho1d[2][k]*dz;
    }
  }
}

*/

#endif
