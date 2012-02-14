/* -*- c++ -*- ----------------------------------------------------------
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

#include "modify.h"
#include "error.h"

namespace GPU_EXTRA {

  inline void check_flag(int error_flag, LAMMPS_NS::Error *error,
                         MPI_Comm &world) { 
    int all_success;
    MPI_Allreduce(&error_flag, &all_success, 1, MPI_INT, MPI_MIN, world);
    if (all_success != 0) {
      if (all_success == -1)
	error->all(FLERR,"Accelerated style in input script but no fix gpu"); 
      else if (all_success == -2)
	error->all(FLERR,
		   "Could not find/initialize a specified accelerator device");
      else if (all_success == -3)
	error->all(FLERR,"Insufficient memory on accelerator");
      else if (all_success == -4)
	error->all(FLERR,"GPU library not compiled for this accelerator");
      else if (all_success == -5)
	error->all(FLERR,
		   "Double precision is not supported on this accelerator");
      else if (all_success == -6)
	error->all(FLERR,"Unable to initialize accelerator for use");
      else if (all_success == -7)
	error->all(FLERR,
                   "Accelerator sharing is not currently supported on system");
      else if (all_success == -8)
	error->all(FLERR,
                   "GPU particle split must be set to 1 for this pair style.");
      else
	error->all(FLERR,"Unknown error in GPU library");
    }
  };

  inline void gpu_ready(LAMMPS_NS::Modify *modify, LAMMPS_NS::Error *error) {
    int ifix = modify->find_fix("package_gpu");
    if (ifix < 0)
      error->all(FLERR,"The package gpu command is required for gpu styles");
  };
}

#endif

/* ERROR/WARNING messages:

E: Accelerated style in input script but no fix gpu

UNDOCUMENTED

E: Could not find/initialize a specified accelerator device

UNDOCUMENTED

E: Insufficient memory on accelerator

UNDOCUMENTED

E: GPU library not compiled for this accelerator

UNDOCUMENTED

E: Double precision is not supported on this accelerator

UNDOCUMENTED

E: Unable to initialize accelerator for use

UNDOCUMENTED

E: Accelerator sharing is not currently supported on system

UNDOCUMENTED

E: GPU particle split must be set to 1 for this pair style.

UNDOCUMENTED

E: Unknown error in GPU library

UNDOCUMENTED

E: The package gpu command is required for gpu styles

UNDOCUMENTED

*/
