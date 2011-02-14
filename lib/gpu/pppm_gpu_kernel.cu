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

#define MAX_STENCIL 8
#define BLOCK_1D 64

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

__device__ inline void atomicFloatAdd(double* address, double val) {
  double old = *address, assumed;
  do { 
    assumed = old;
    old = __longlong_as_double( atomicCAS((unsigned long long int*)address, 
                                          __double_as_longlong(assumed),
                                          __double_as_longlong(val +
                                          assumed)));
  } while (assumed != old); 
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

__device__ inline void atomicFloatAdd(float *address, float val)
{
       int i_val = __float_as_int(val);
       int tmp0 = 0;
       int tmp1;

       while( (tmp1 = atomicCAS((int *)address, tmp0, i_val)) != tmp0)
       {
               tmp0 = tmp1;
               i_val = __float_as_int(val + __int_as_float(tmp1));
       }
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
                           const numtyp b_lo_x, const numtyp b_lo_y,
                           const numtyp b_lo_z, const numtyp delxinv,
                           const numtyp delyinv, const numtyp delzinv,
                           const int nlocal_x, const int nlocal_y,
                           const int nlocal_z, const int atom_stride,
                           const int max_atoms, __global int *error) {
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;
  int nx,ny,nz;
  numtyp tx,ty,tz;

  if (ii<nlocal) {
    numtyp4 p=fetch_pos(ii,x_);

    tx=(p.x-b_lo_x)*delxinv;
    nx=int(tx);
    ty=(p.y-b_lo_y)*delyinv;
    ny=int(ty);
    tz=(p.z-b_lo_z)*delzinv;
    nz=int(tz);

    if (tx<0 || ty<0 || tz<0 || nx>=nlocal_x || ny>=nlocal_y || nz>=nlocal_z)
      *error=1;
    else {
      int i=nz*nlocal_y*nlocal_x+ny*nlocal_x+nx;
      int old=atom_add(counts+i, 1);
      if (old==max_atoms) {
        *error=2;
        atom_add(counts+i,-1);
      } else
        ans[atom_stride*old+i]=ii;
    }
  }
}

__kernel void make_rho(__global numtyp4 *x_, __global numtyp *q_,
                       __global int *counts, __global int *atoms,
                       __global numtyp *brick, __global numtyp *_rho_coeff,
                       const int atom_stride, const int npts_x,
                       const int npts_y, const int nlocal_x, const int nlocal_y,
                       const int nlocal_z, const numtyp b_lo_x,
                       const numtyp b_lo_y, const numtyp b_lo_z,
                       const numtyp delxinv, const numtyp delyinv,
                       const numtyp delzinv, const int order,
                       const numtyp delvolinv) {
  __local numtyp rho_coeff[MAX_STENCIL*MAX_STENCIL];
  int nx=THREAD_ID_X;
  int ny=THREAD_ID_Y;
  if (nx<order && ny<order) {
    int ri=nx*order+ny;
    rho_coeff[ri]=_rho_coeff[ri];
  }
  __syncthreads();
  
  nx+=BLOCK_ID_X*BLOCK_SIZE_X;
  ny+=BLOCK_ID_Y*BLOCK_SIZE_Y;
  int nz=0;
  
  if (nx<nlocal_x && ny<nlocal_y) {
    int z_stride=nlocal_x*nlocal_y;
    int z_pos=nz*z_stride+ny*nlocal_x+nx;
    for ( ; nz<nlocal_z; nz++) {
      int natoms=counts[z_pos];
      for (int row=0; row<natoms; row++) {
        int atom=atoms[atom_stride*row+z_pos];
        numtyp4 p=fetch_pos(atom,x_);
        numtyp z0=delvolinv*fetch_q(atom,q_);
        
        numtyp dx = nx - (p.x-b_lo_x)*delxinv;
        numtyp dy = ny - (p.y-b_lo_y)*delyinv;
        numtyp dz = nz - (p.z-b_lo_z)*delzinv;

        numtyp rho1d[2][MAX_STENCIL];
        for (int k = 0; k < order; k++) {
          rho1d[0][k] = (numtyp)0.0;
          rho1d[1][k] = (numtyp)0.0;
          for (int l = order-1; l >= 0; l--) {
            rho1d[0][k] = rho_coeff[l*order+k] + rho1d[0][k]*dx;
            rho1d[1][k] = rho_coeff[l*order+k] + rho1d[1][k]*dy;
          }
        }
        
        for (int n = 0; n < order; n++) {
          numtyp rho1d_2 = (numtyp)0.0;
          for (int k = order-1; k >= 0; k--)
            rho1d_2 = rho_coeff[k*order+n] + rho1d_2*dz;
          numtyp y0 = z0*rho1d_2;
          int mz = (n+nz)*npts_y*npts_x + ny*npts_x +nx;
          for (int m = 0; m < order; m++) {
	          numtyp x0 = y0*rho1d[1][m];
	          for (int l = 0; l < order; l++) {
              atomicFloatAdd(brick+mz+l,x0*rho1d[0][l]);
	          }
	          mz+=npts_x;
	        }
	      }
	    }
	    z_pos+=z_stride;
	  }
	}
}

/* --------------------------- */

__kernel void make_rho2(__global numtyp4 *x_, __global numtyp *q_,
                       __global int *counts, __global int *atoms,
                       __global numtyp *brick, __global numtyp *_rho_coeff,
                       const int atom_stride, const int npts_x,
                       const int npts_y, const int npts_z, 
                       const int nlocal_x, const int nlocal_y,
                       const int nlocal_z, const int nlow2,
                       const numtyp b_lo_x,
                       const numtyp b_lo_y, const numtyp b_lo_z,
                       const numtyp delxinv, const numtyp delyinv,
                       const numtyp delzinv, const int order,
                       const numtyp delvolinv) {
  __local numtyp rho_coeff[MAX_STENCIL*MAX_STENCIL];
  __local numtyp front[BLOCK_1D+MAX_STENCIL];
  __local int nx,ny,x_start,y_start,x_stop,y_stop;
  __local int z_stride, z_local_stride;

  int tx=THREAD_ID_X;
  int tx_halo=BLOCK_1D+tx;
  if (tx<order*order)
    rho_coeff[tx]=_rho_coeff[tx];
    
  if (tx==0) {
    nx=BLOCK_ID_X;
    ny=BLOCK_ID_Y;
    x_start=0;
    y_start=0;
    x_stop=order;
    y_stop=order;
    if (nx<nlow2)
      x_start=nlow2-nx;
    if (ny<nlow2)
      y_start=nlow2-ny;
    if (nx>=nlocal_x)
      x_stop-=nx-nlocal_x+1;
    if (ny>=nlocal_y)
      y_stop-=ny-nlocal_y+1;
    z_stride=npts_x*npts_y*BLOCK_1D;
    z_local_stride=nlocal_x*nlocal_y*BLOCK_1D;
  }
  
  if (tx<order) 
    front[tx_halo]=(numtyp)0.0;
    
  __syncthreads();

  numtyp ans[MAX_STENCIL];
  int loop_count=npts_z/BLOCK_1D+1;
  int nz=tx;
  int pt = nz*npts_x*npts_y + ny*npts_x + nx;
  int z_local = nz*nlocal_x*nlocal_y;
  for (int i=0 ; i<loop_count; i++) {
    int offset1,offset2;
    if (i%2) {
      offset1=0;
      offset2=1;
    } else {
      offset1=1;
      offset2=0;
    }
    for (int n=0; n<order; n++)
      ans[n]=(numtyp)0.0;
    if (nz<nlocal_z) {
      for (int m=y_start; m<y_stop; m++) {
        int y_pos=ny+m-nlow2;
        int y_local=y_pos*nlocal_x;
        for (int l=x_start; l<x_stop; l++) {
          int x_pos=nx+l-nlow2;
          int pos=z_local+y_local+x_pos;
          int natoms=counts[pos]*atom_stride;
          for (int row=pos; row<natoms; row+=atom_stride) {
            int atom=atoms[row];
            numtyp4 p=fetch_pos(atom,x_);
            numtyp z0=delvolinv*fetch_q(atom,q_);
      
            numtyp dx = x_pos - (p.x-b_lo_x)*delxinv;
            numtyp dy = y_pos - (p.y-b_lo_y)*delyinv;
            numtyp dz = nz - (p.z-b_lo_z)*delzinv;
            
            numtyp rho1d_1 = (numtyp)0.0;
            numtyp rho1d_0 = (numtyp)0.0;
            for (int k = order; k > 0; k--) {
              rho1d_1 = rho_coeff[k*order-m-1] + rho1d_1*dy;
              rho1d_0 = rho_coeff[k*order-l-1] + rho1d_0*dx;
            }
            z0*=rho1d_1*rho1d_0;

            for (int n=0; n<order; n++) {
              numtyp rho1d_2 = (numtyp)0.0;
              for (int k = order-1; k >= 0; k--)
                rho1d_2 = rho_coeff[k*order+n] + rho1d_2*dz;
              ans[n]+=z0*rho1d_2;
            }
          }
        }
      }
    }
    
    __syncthreads();
    if (tx<order) {
      front[tx]=front[tx_halo];
      front[tx_halo]=(numtyp)0.0;
    } else 
      front[tx]=(numtyp)0.0;
    
    for (int n=0; n<order; n++) {
      front[tx+n]+=ans[n];
      __syncthreads();
    }

    if (nz<npts_z)
      brick[pt]=front[tx];
    nz+=BLOCK_1D;
    pt+=z_stride;
    z_local+=z_local_stride;
  }
}

/* --------------------------- */

__kernel void make_rho3(__global numtyp4 *x_, __global numtyp *q_,
                       const int nlocal,
                       __global numtyp *brick, __global numtyp *_rho_coeff,
                       const int npts_x,
                       const int npts_y, const int nlocal_x, const int nlocal_y,
                       const int nlocal_z, const numtyp b_lo_x,
                       const numtyp b_lo_y, const numtyp b_lo_z,
                       const numtyp delxinv, const numtyp delyinv,
                       const numtyp delzinv, const numtyp shift,
                       const int order,
                       const numtyp delvolinv, __global int *error) {
  __local numtyp rho_coeff[MAX_STENCIL*MAX_STENCIL];
  int ii=THREAD_ID_X;
  if (ii<order*order)
    rho_coeff[ii]=_rho_coeff[ii];
  __syncthreads();
  
  ii+=BLOCK_ID_X*BLOCK_SIZE_X;
  
  int nx,ny,nz;
  numtyp tx,ty,tz;

  if (ii<nlocal) {
    numtyp4 p=fetch_pos(ii,x_);

    tx=(p.x-b_lo_x)*delxinv;
    nx=int(tx);
    ty=(p.y-b_lo_y)*delyinv;
    ny=int(ty);
    tz=(p.z-b_lo_z)*delzinv;
    nz=int(tz);

    if (tx<0 || ty<0 || tz<0 || nx>=nlocal_x || ny>=nlocal_y || nz>=nlocal_z)
      *error=1;
    else {
      numtyp z0=delvolinv*fetch_q(ii,q_);
        
      numtyp dx = nx+shift - tx;
      numtyp dy = ny+shift - ty;
      numtyp dz = nz+shift - tz;

      numtyp rho1d[2][MAX_STENCIL];
      for (int k = 0; k < order; k++) {
        rho1d[0][k] = (numtyp)0.0;
        rho1d[1][k] = (numtyp)0.0;
        for (int l = order-1; l >= 0; l--) {
          rho1d[0][k] = rho_coeff[l*order+k] + rho1d[0][k]*dx;
          rho1d[1][k] = rho_coeff[l*order+k] + rho1d[1][k]*dy;
        }
      }
        
      for (int n = 0; n < order; n++) {
        numtyp rho1d_2 = (numtyp)0.0;
        for (int k = order-1; k >= 0; k--)
          rho1d_2 = rho_coeff[k*order+n] + rho1d_2*dz;
        numtyp y0 = z0*rho1d_2;
        int mz = (n+nz)*npts_y*npts_x + ny*npts_x +nx;
        for (int m = 0; m < order; m++) {
          numtyp x0 = y0*rho1d[1][m];
	        for (int l = 0; l < order; l++) {
            atomicFloatAdd(brick+mz+l,x0*rho1d[0][l]);
	        }
          mz+=npts_x;
        }
	    }
	  }
	}
}

#endif

