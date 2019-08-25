// **************************************************************************
//                                  atom.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for atom data casting
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
#endif

__kernel void kernel_cast_x(__global numtyp4 *restrict x_type,
                            const __global double *restrict x,
                            const __global int *restrict type,
                            const int nall) {
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
