// **************************************************************************
//                                  pppm.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for PPPM acceleration
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : brownw@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL

#include "lal_preprocessor.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
texture<float> q_tex;
#else
texture<int4,1> pos_tex;
texture<int2> q_tex;
#endif

// Allow PPPM to compile without atomics for NVIDIA 1.0 cards, error
// generated at runtime with use of pppm/gpu
#if (__CUDA_ARCH__ < 110)
#define atomicAdd(x,y) *(x)+=0
#endif

#else
#define pos_tex x_
#define q_tex q_
#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics: enable
#endif

// Number of threads per pencil for charge spread
#define PENCIL_SIZE MEM_THREADS
// Number of pencils per block for charge spread
#define BLOCK_PENCILS (PPPM_BLOCK_1D/PENCIL_SIZE)

__kernel void particle_map(const __global numtyp4 *restrict x_,  
                           const __global numtyp *restrict q_,
                           const grdtyp delvolinv, const int nlocal, 
                           __global int *restrict counts, 
                           __global grdtyp4 *restrict ans, 
                           const grdtyp b_lo_x, const grdtyp b_lo_y,
                           const grdtyp b_lo_z, const grdtyp delxinv,
                           const grdtyp delyinv, const grdtyp delzinv,
                           const int nlocal_x, const int nlocal_y,
                           const int nlocal_z, const int atom_stride,
                           const int max_atoms, 
                           __global int *restrict error) {
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;

  // Resequence the atom indices to avoid collisions during atomic ops
  int nthreads=GLOBAL_SIZE_X;
  ii=fast_mul(ii,PPPM_BLOCK_1D);
  ii-=(ii/nthreads)*(nthreads-1);

  int nx,ny,nz;

  if (ii<nlocal) {
    numtyp4 p;
    fetch4(p,ii,pos_tex);
    grdtyp4 delta;
    fetch(delta.w,ii,q_tex);
    delta.w*=delvolinv;
    
    if (delta.w!=(grdtyp)0.0) {
      delta.x=(p.x-b_lo_x)*delxinv;
      nx=delta.x;
      delta.y=(p.y-b_lo_y)*delyinv;
      ny=delta.y;
      delta.z=(p.z-b_lo_z)*delzinv;
      nz=delta.z;

      if (delta.x<(grdtyp)0 || delta.y<(grdtyp)0 || delta.z<(grdtyp)0 || 
          nx>=nlocal_x || ny>=nlocal_y || nz>=nlocal_z)
        *error=1;
      else {
        delta.x=nx+(grdtyp)0.5-delta.x;
        delta.y=ny+(grdtyp)0.5-delta.y;
        delta.z=nz+(grdtyp)0.5-delta.z;
      
        int i=nz*nlocal_y*nlocal_x+ny*nlocal_x+nx;
        int old=atom_add(counts+i, 1);
        if (old>=max_atoms) {
          *error=2;
          atom_add(counts+i, -1);
        } else
          ans[atom_stride*old+i]=delta;
      }
    }
  }
}

/* --------------------------- */

__kernel void make_rho(const __global int *restrict counts, 
                       const __global grdtyp4 *restrict atoms,
                       __global grdtyp *restrict brick, 
                       const __global grdtyp *restrict _rho_coeff,
                       const int atom_stride, const int npts_x,
                       const int npts_y, const int npts_z, const int nlocal_x,
                       const int nlocal_y, const int nlocal_z,
                       const int order_m_1, const int order, const int order2) {
  __local grdtyp rho_coeff[PPPM_MAX_SPLINE*PPPM_MAX_SPLINE];
  __local grdtyp front[BLOCK_PENCILS][PENCIL_SIZE+PPPM_MAX_SPLINE];
  __local grdtyp ans[PPPM_MAX_SPLINE][PPPM_BLOCK_1D];
  
  int tid=THREAD_ID_X;
  if (tid<order2+order)
    rho_coeff[tid]=_rho_coeff[tid];
    
  int pid=tid/PENCIL_SIZE;
  int fid=tid%PENCIL_SIZE;
  int fid_halo=PENCIL_SIZE+fid;
  if (fid<order) 
    front[pid][fid_halo]=(grdtyp)0.0;

  __syncthreads();

  int bt=BLOCK_ID_X*BLOCK_PENCILS+pid;
  int ny=bt%npts_y;
  int nz=bt/npts_y;
  int y_start=0;
  int z_start=0;
  int y_stop=order;
  int z_stop=order;
  if (ny<order_m_1)
    y_start=order_m_1-ny;
  if (nz<order_m_1)
    z_start=order_m_1-nz;
  if (ny>=nlocal_y)
    y_stop-=ny-nlocal_y+1;
  if (nz>=nlocal_z)
    z_stop-=nz-nlocal_z+1;
  int z_stride=fast_mul(nlocal_x,nlocal_y);

  int loop_count=npts_x/PENCIL_SIZE+1;
  int nx=fid;
  int pt=fast_mul(nz,fast_mul(npts_y,npts_x))+fast_mul(ny,npts_x)+nx;
  for (int i=0 ; i<loop_count; i++) {
    for (int n=0; n<order; n++)
      ans[n][tid]=(grdtyp)0.0;
    if (nx<nlocal_x && nz<npts_z) {
      int z_pos=fast_mul(nz+z_start-order_m_1,z_stride);
      for (int m=z_start; m<z_stop; m++) {
        int y_pos=fast_mul(ny+y_start-order_m_1,nlocal_x);
        for (int l=y_start; l<y_stop; l++) {
          int pos=z_pos+y_pos+nx;
          int natoms=fast_mul(counts[pos],atom_stride);
          for (int row=pos; row<natoms; row+=atom_stride) {
            grdtyp4 delta=atoms[row];
      
            grdtyp rho1d_1=(grdtyp)0.0;
            grdtyp rho1d_2=(grdtyp)0.0;
            for (int k=order2+order-1; k > -1; k-=order) {
              rho1d_1=rho_coeff[k-l]+rho1d_1*delta.y;
              rho1d_2=rho_coeff[k-m]+rho1d_2*delta.z;
            }
            delta.w*=rho1d_1*rho1d_2;

            for (int n=0; n<order; n++) {
              grdtyp rho1d_0=(grdtyp)0.0;
              for (int k=order2+n; k>=n; k-=order)
                rho1d_0=rho_coeff[k]+rho1d_0*delta.x;
              ans[n][tid]+=delta.w*rho1d_0;
            }
          }
          y_pos+=nlocal_x;
        }
        z_pos+=z_stride;
      }
    }
    
    __syncthreads();
    if (fid<order) {
      front[pid][fid]=front[pid][fid_halo];
      front[pid][fid_halo]=(grdtyp)0.0;
    } else 
      front[pid][fid]=(grdtyp)0.0;
    
    for (int n=0; n<order; n++) {
      front[pid][fid+n]+=ans[n][tid];
      __syncthreads();
    }

    if (nx<npts_x && nz<npts_z)
      brick[pt]=front[pid][fid];
    pt+=PENCIL_SIZE;
    nx+=PENCIL_SIZE;
  }
}

__kernel void interp(const __global numtyp4 *restrict x_, 
                     const __global numtyp *restrict q_,
                     const int nlocal, 
                     const __global grdtyp4 *restrict brick,
                     const __global grdtyp *restrict _rho_coeff, 
                     const int npts_x, const int npts_yx, const grdtyp b_lo_x,
                     const grdtyp b_lo_y, const grdtyp b_lo_z,
                     const grdtyp delxinv,  const grdtyp delyinv,
                     const grdtyp delzinv, const int order,
                     const int order2, const grdtyp qqrd2e_scale, 
                     __global acctyp4 *restrict ans) {
  __local grdtyp rho_coeff[PPPM_MAX_SPLINE*PPPM_MAX_SPLINE];
  __local grdtyp rho1d_0[PPPM_MAX_SPLINE][PPPM_BLOCK_1D];
  __local grdtyp rho1d_1[PPPM_MAX_SPLINE][PPPM_BLOCK_1D];

  int tid=THREAD_ID_X;
  if (tid<order2+order)
    rho_coeff[tid]=_rho_coeff[tid];
  __syncthreads();
  
  int ii=tid+BLOCK_ID_X*BLOCK_SIZE_X;
  
  int nx,ny,nz;
  grdtyp tx,ty,tz;

  if (ii<nlocal) {
    numtyp4 p;
    fetch4(p,ii,pos_tex);
    grdtyp qs;
    fetch(qs,ii,q_tex);
    qs*=qqrd2e_scale;

    acctyp4 ek;
    ek.x=(acctyp)0.0;
    ek.y=(acctyp)0.0;
    ek.z=(acctyp)0.0;
    if (qs!=(grdtyp)0.0) {
      tx=(p.x-b_lo_x)*delxinv;
      nx=tx;
      ty=(p.y-b_lo_y)*delyinv;
      ny=ty;
      tz=(p.z-b_lo_z)*delzinv;
      nz=tz;

      grdtyp dx=nx+(grdtyp)0.5-tx;
      grdtyp dy=ny+(grdtyp)0.5-ty;
      grdtyp dz=nz+(grdtyp)0.5-tz;

      for (int k=0; k<order; k++) {
        rho1d_0[k][tid]=(grdtyp)0.0;
        rho1d_1[k][tid]=(grdtyp)0.0;
        for (int l=order2+k; l>=k; l-=order) {
          rho1d_0[k][tid]=rho_coeff[l]+rho1d_0[k][tid]*dx;
          rho1d_1[k][tid]=rho_coeff[l]+rho1d_1[k][tid]*dy;
        }
      }
        
      int mz=fast_mul(nz,npts_yx)+nx;
      for (int n=0; n<order; n++) {
        grdtyp rho1d_2=(grdtyp)0.0;
        for (int k=order2+n; k>=n; k-=order)
          rho1d_2=rho_coeff[k]+rho1d_2*dz;
        grdtyp z0=qs*rho1d_2;
        int my=mz+fast_mul(ny,npts_x);
        for (int m=0; m<order; m++) {
          grdtyp y0=z0*rho1d_1[m][tid];
  	      for (int l=0; l<order; l++) {
  	        grdtyp x0=y0*rho1d_0[l][tid];
  	        grdtyp4 el=brick[my+l];
  	        ek.x-=x0*el.x;
  	        ek.y-=x0*el.y;
  	        ek.z-=x0*el.z;
  	      }
          my+=npts_x;
        }
        mz+=npts_yx;
  	  }
    }
    ans[ii]=ek;
	}
}

