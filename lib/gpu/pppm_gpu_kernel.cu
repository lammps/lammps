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
#define GLOBAL_SIZE_X get_global_size(0)
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

  // Resequence the atom indices to avoid collisions during atomic ops
  int nthreads=GLOBAL_SIZE_X;
  ii=mul24(ii,BLOCK_1D);
  ii-=int(ii/nthreads)*(nthreads-1);

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

/* --------------------------- */

__kernel void make_rho(__global numtyp4 *x_, __global numtyp *q_,
                       __global int *counts, __global int *atoms,
                       __global numtyp *brick, __global numtyp *_rho_coeff,
                       const int atom_stride, const int npts_x,
                       const int npts_yx, const int npts_z, const int nlocal_x,
                       const int nlocal_y, const int nlocal_z,
                       const int order_m_1, const numtyp b_lo_x,
                       const numtyp b_lo_y, const numtyp b_lo_z,
                       const numtyp delxinv, const numtyp delyinv,
                       const numtyp delzinv, const int order, const int order2,
                       const numtyp delvolinv) {
  __local numtyp rho_coeff[MAX_STENCIL*MAX_STENCIL];
  __local numtyp front[BLOCK_1D+MAX_STENCIL];
  __local numtyp ans[MAX_STENCIL][BLOCK_1D];
  __local int nx,ny,x_start,y_start,x_stop,y_stop;
  __local int z_stride, z_local_stride;

  int tx=THREAD_ID_X;
  int tx_halo=BLOCK_1D+tx;
  if (tx<order2+order)
    rho_coeff[tx]=_rho_coeff[tx];
    
  if (tx==0) {
    nx=BLOCK_ID_X;
    ny=BLOCK_ID_Y;
    x_start=0;
    y_start=0;
    x_stop=order;
    y_stop=order;
    if (nx<order_m_1)
      x_start=order_m_1-nx;
    if (ny<order_m_1)
      y_start=order_m_1-ny;
    if (nx>=nlocal_x)
      x_stop-=nx-nlocal_x+1;
    if (ny>=nlocal_y)
      y_stop-=ny-nlocal_y+1;
    z_stride=mul24(npts_yx,BLOCK_1D);
    z_local_stride=mul24(mul24(nlocal_x,nlocal_y),BLOCK_1D);
  }
  
  if (tx<order) 
    front[tx_halo]=(numtyp)0.0;
    
  __syncthreads();

  int loop_count=npts_z/BLOCK_1D+1;
  int nz=tx;
  int pt=mul24(nz,npts_yx)+mul24(ny,npts_x)+nx;
  int z_local=mul24(mul24(nz,nlocal_x),nlocal_y);
  for (int i=0 ; i<loop_count; i++) {
    for (int n=0; n<order; n++)
      ans[n][tx]=(numtyp)0.0;
    if (nz<nlocal_z) {
      for (int m=y_start; m<y_stop; m++) {
        int y_pos=ny+m-order_m_1;
        int y_local=mul24(y_pos,nlocal_x);
        for (int l=x_start; l<x_stop; l++) {
          int x_pos=nx+l-order_m_1;
          int pos=z_local+y_local+x_pos;
          int natoms=mul24(counts[pos],atom_stride);
          for (int row=pos; row<natoms; row+=atom_stride) {
            int atom=atoms[row];
            numtyp4 p=fetch_pos(atom,x_);
            numtyp z0=delvolinv*fetch_q(atom,q_);
      
            numtyp dx=x_pos-(p.x-b_lo_x)*delxinv;
            numtyp dy=y_pos-(p.y-b_lo_y)*delyinv;
            numtyp dz=nz-(p.z-b_lo_z)*delzinv;
            
            numtyp rho1d_1=(numtyp)0.0;
            numtyp rho1d_0=(numtyp)0.0;
            for (int k=order2+order-1; k > -1; k-=order) {
              rho1d_1=rho_coeff[k-m]+rho1d_1*dy;
              rho1d_0=rho_coeff[k-l]+rho1d_0*dx;
            }
            z0*=rho1d_1*rho1d_0;

            for (int n=0; n<order; n++) {
              numtyp rho1d_2=(numtyp)0.0;
              for (int k=order2+n; k>=n; k-=order)
                rho1d_2=rho_coeff[k]+rho1d_2*dz;
              ans[n][tx]+=z0*rho1d_2;
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
      front[tx+n]+=ans[n][tx];
      __syncthreads();
    }

    if (nz<npts_z)
      brick[pt]=front[tx];
    nz+=BLOCK_1D;
    pt+=z_stride;
    z_local+=z_local_stride;
  }
}

__kernel void field_force(__global numtyp4 *x_, __global numtyp *q_,
                          const int nlocal, __global numtyp *x_brick,
                          __global numtyp *y_brick, __global numtyp *z_brick,
                          __global numtyp *_rho_coeff, const int npts_x,
                          const int npts_yx, const numtyp b_lo_x,
                          const numtyp b_lo_y, const numtyp b_lo_z,
                          const numtyp delxinv,  const numtyp delyinv,
                          const numtyp delzinv, const int order,
                          const int order2, const numtyp qqrd2e, 
                          const numtyp scale, __global acctyp4 *ans) {
  __local numtyp rho_coeff[MAX_STENCIL*MAX_STENCIL];
  __local numtyp rho1d_0[MAX_STENCIL][BLOCK_1D];
  __local numtyp rho1d_1[MAX_STENCIL][BLOCK_1D];

  int tid=THREAD_ID_X;
  if (tid<order2+order)
    rho_coeff[tid]=_rho_coeff[tid];
  __syncthreads();
  
  int ii=tid+BLOCK_ID_X*BLOCK_SIZE_X;
  
  int nx,ny,nz;
  numtyp tx,ty,tz;

  if (ii<nlocal) {
    numtyp4 p=fetch_pos(ii,x_);
    numtyp qs=qqrd2e*scale*fetch_q(ii,q_);

    tx=(p.x-b_lo_x)*delxinv;
    nx=int(tx);
    ty=(p.y-b_lo_y)*delyinv;
    ny=int(ty);
    tz=(p.z-b_lo_z)*delzinv;
    nz=int(tz);

    numtyp dx=nx+(numtyp)0.5-tx;
    numtyp dy=ny+(numtyp)0.5-ty;
    numtyp dz=nz+(numtyp)0.5-tz;

    for (int k=0; k<order; k++) {
      rho1d_0[k][tid]=(numtyp)0.0;
      rho1d_1[k][tid]=(numtyp)0.0;
      for (int l=order2+k; l>=k; l-=order) {
        rho1d_0[k][tid]=rho_coeff[l]+rho1d_0[k][tid]*dx;
        rho1d_1[k][tid]=rho_coeff[l]+rho1d_1[k][tid]*dy;
      }
    }
        
    numtyp4 ek;
    ek.x=(acctyp)0.0;
    ek.y=(acctyp)0.0;
    ek.z=(acctyp)0.0;
    int mz=mul24(nz,npts_yx)+nx;
    for (int n=0; n<order; n++) {
      numtyp rho1d_2=(numtyp)0.0;
      for (int k=order2+n; k>=n; k-=order)
        rho1d_2=rho_coeff[k]+rho1d_2*dz;
      numtyp z0=qs*rho1d_2;
      int my=mz+mul24(ny,npts_x);
      for (int m=0; m<order; m++) {
        numtyp y0=z0*rho1d_1[m][tid];
	      for (int l=0; l<order; l++) {
	        numtyp x0=y0*rho1d_0[l][tid];
	        ek.x-=x0*x_brick[my+l];
	        ek.y-=x0*y_brick[my+l];
	        ek.z-=x0*z_brick[my+l];
	      }
        my+=npts_x;
      }
      mz+=npts_yx;
	  }
    ans[ii]=ek;
	}
}

#endif

