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

#ifdef NV_KERNEL

#include "geryon/ucl_nv_kernel.h"

#else

#define GLOBAL_ID_X get_global_id(0)

#endif

#define O2_BLOCK_SIZE 64

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp4 double4
#else
#define numtyp float
#define numtyp4 float4
#endif

__kernel void kernel_unpack(__global int *dev_nbor, __global int *dev_ij,
                            const int inum) {
  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;

  if (ii<inum) {
    __global int *nbor=dev_nbor+ii+inum;
    int numj=*nbor;
    nbor+=inum;
    __global int *list=dev_ij+*nbor;
    __global int *list_end=list+numj;
  
    for ( ; list<list_end; list++) {
      *nbor=*list;
      nbor+=inum;
    }
  } // if ii
}

__kernel void kernel_nbor_o2(__global numtyp4 *x_, __global int *dev_nbor, 
                             __global int *host_nbor_list, const int inum,
                             const int nt, const int nall, 
                             numtyp cell_size_sq) {
  int ii=THREAD_ID_X;
  int i=ii+mul24((int)BLOCK_ID_X,(int)BLOCK_SIZE_X);
  int numj=0;

  int *nbor;
  int nbor_pitch;
  
  if (i < inum) {
    nbor_pitch=inum;
    nbor=dev_nbor+i;
    *nbor=i;
    nbor+=nbor_pitch;
  } else {
    nbor_pitch=nt-inum;
    nbor=host_nbor_list+i-inum;
  }
  nbor+=nbor_pitch;
  
  // Shared memory for atom positions
  __local numtyp4 x[O2_BLOCK_SIZE];

  numtyp4 x_me;
  if (i<nt)
    x_me=x_[i];

  for (int j=0; j<nall; j+=O2_BLOCK_SIZE) {
    int thread_j=j+ii;
    if (thread_j<nall)
      x[ii]=x_[thread_j];

    __syncthreads();
    
    if (i<nt) {
      for (int k=0; k<O2_BLOCK_SIZE; k++) {
        numtyp rsq=x_me.x-x[k].x;
        rsq*=rsq;
        numtyp tsq=x_me.y-x[k].y;
        rsq+=tsq*tsq;
        tsq=x_me.z-x[k].z;
        rsq+=tsq*tsq;
        
        if (rsq<cell_size_sq) {
          int tj=j+k;
          if (tj<nall && tj!=i) {
            if (i<inum || j<inum || j>i) {
              *nbor=j;
              nbor+=nbor_pitch;
              numj++;
            }
          }
        }
      }
    }
  }

  if (i < inum)
    nbor=dev_nbor+i+nbor_pitch;
  else
    nbor=host_nbor_list+i-inum;
  if (i < nt)
    *nbor=numj;
}
