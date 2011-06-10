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

