// **************************************************************************
//                                  atom.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for handling CPU generated neighbor lists
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : brownw@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "preprocessor.h"
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

