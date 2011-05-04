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
   Contributing author: Mike Brown (ORNL)
------------------------------------------------------------------------- */

#ifndef LMP_GPU_EXTRA_H
#define LMP_GPU_EXTRA_H

#include "error.h"

namespace GPU_EXTRA {

  inline void check_flag(int error_flag, LAMMPS_NS::Error *error,
                         MPI_Comm &world) { 
    int all_success;
    MPI_Allreduce(&error_flag, &all_success, 1, MPI_INT, MPI_MIN, world);
    if (all_success != 0) {
      if (all_success == -1)
	error->all("Accelerated style in input script but no fix gpu."); 
      else if (all_success == -2)
	error->all("Could not find/initialize a specified accelerator device.");
      else if (all_success == -3)
	error->all("Insufficient memory on accelerator.");
      else if (all_success == -4)
	error->all("GPU library not compiled for this accelerator.");
      else if (all_success == -5)
	error->all("Double precision is not supported on this accelerator.");
      else
	error->all("Unknown error in GPU library.");
    }
  }

}

#endif
