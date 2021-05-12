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

#ifndef LMP_MDI_INTERFACE_H
#define LMP_MDI_INTERFACE_H

#if defined(LMP_USER_MDI)

// true interface to MDI
// used when MDI is installed

#include "mdi.h"

#else

// dummy interface to MDI
// needed for compiling when MDI is not installed

  typedef int MDI_Comm;
  int MDI_Init(int* argc, char ***argv) {return 0;};
  int MDI_Initialized(int* flag) {return 0;};
  int MDI_MPI_get_world_comm(void* world_comm) {return 0;};
  int MDI_Plugin_get_argc(int* argc) {return 0;};
  int MDI_Plugin_get_argv(char*** argv) {return 0;};

#endif
#endif
