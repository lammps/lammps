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

#ifndef PAIR_GPU_KERNEL_H
#define PAIR_GPU_KERNEL_H

#ifdef NV_KERNEL

#include "nv_kernel_def.h"

#else

#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define GLOBAL_ID_X get_global_id(0)
#define THREAD_ID_X get_local_id(0)
#define BLOCK_ID_X get_group_id(0)
#define BLOCK_SIZE_X get_local_size(0)
#define __syncthreads() barrier(CLK_LOCAL_MEM_FENCE)
#define MAX_SHARED_TYPES 8

#endif

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp2 double2
#define numtyp4 double4
#else
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#endif

// ---------------------------------------------------------------------------
// Unpack neighbors from dev_ij array into dev_nbor matrix for coalesced access
// -- Only unpack neighbors matching the specified inclusive range of forms
// -- Only unpack neighbors within cutoff
// ---------------------------------------------------------------------------
__kernel void kernel_gb_nbor(__global numtyp4 *x_, __global numtyp2 *cut_form, 
                             const int ntypes, __global int *dev_nbor,
                             const int nbor_pitch, 
                             const int start, const int inum, 
                             __global int *dev_ij, const int form_low, 
                             const int form_high, const int nall) {
                                
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X+start;

  if (ii<inum) {
    __global int *nbor=dev_ij+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    __global int *list_end=nbor+mul24(numj,nbor_pitch);
    __global int *packed=dev_nbor+ii+nbor_pitch+nbor_pitch;
  
    numtyp4 ix=x_[i];
    int iw=ix.w;
    int itype=mul24(iw,ntypes);
    int newj=0;  
    for ( ; nbor<list_end; nbor+=nbor_pitch) {
      int j=*nbor;
      if (j>=nall)
        j%=nall;
      numtyp4 jx=x_[j];
      int jtype=jx.w;
      int mtype=itype+jtype;
      numtyp2 cf=cut_form[mtype];
      if (cf.y>=form_low && cf.y<=form_high) {
        // Compute r12;
        numtyp rsq=jx.x-ix.x;
        rsq*=rsq;
        numtyp t=jx.y-ix.y;
        rsq+=t*t;
        t=jx.z-ix.z;
        rsq+=t*t;

        if (rsq<cf.x) {
          *packed=j;
          packed+=nbor_pitch;
          newj++;
        }
      }
    }
    dev_nbor[ii+nbor_pitch]=newj;
  }
}

// ---------------------------------------------------------------------------
// Unpack neighbors from dev_ij array into dev_nbor matrix for coalesced access
// -- Only unpack neighbors matching the specified inclusive range of forms
// -- Only unpack neighbors within cutoff
// -- Fast version of routine that uses shared memory for LJ constants
// ---------------------------------------------------------------------------
__kernel void kernel_gb_nbor_fast(__global numtyp4 *x_, 
                                  __global numtyp2 *cut_form,
                                  __global int *dev_nbor, 
                                  const int nbor_pitch, 
                                  const int start, const int inum, 
                                  __global int *dev_ij, const int form_low, 
                                  const int form_high, const int nall) {
                                
  int ii=THREAD_ID_X;
  __local int form[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (ii<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[ii]=cut_form[ii].x;
    form[ii]=cut_form[ii].y;
  }
  ii+=mul24((int)BLOCK_SIZE_X,(int)BLOCK_ID_X)+start;
  __syncthreads();

  if (ii<inum) {
    __global int *nbor=dev_ij+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    __global int *list_end=nbor+mul24(numj,nbor_pitch);
    __global int *packed=dev_nbor+ii+nbor_pitch+nbor_pitch;
  
    numtyp4 ix=x_[i];
    int iw=ix.w;
    int itype=mul24((int)MAX_SHARED_TYPES,iw);

    int newj=0;  
    for ( ; nbor<list_end; nbor+=nbor_pitch) {
      int j=*nbor;
      if (j>=nall)
        j%=nall;
      numtyp4 jx=x_[j];
      int jtype=jx.w;
      int mtype=itype+jtype;
      
      if (form[mtype]>=form_low && form[mtype]<=form_high) {
        // Compute r12;
        numtyp rsq=jx.x-ix.x;
        rsq*=rsq;
        numtyp t=jx.y-ix.y;
        rsq+=t*t;
        t=jx.z-ix.z;
        rsq+=t*t;

        if (rsq<cutsq[mtype]) {
          *packed=j;
          packed+=nbor_pitch;
          newj++;
        }
      }
    }
    dev_nbor[ii+nbor_pitch]=newj;
  }
}

#endif
