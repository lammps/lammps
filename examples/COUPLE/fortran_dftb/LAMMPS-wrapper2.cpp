/* -----------------------------------------------------------------------
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    https://www.lammps.org/
    Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories
 
    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under 
    the GNU General Public License.
 
    See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
    Contributing author:  Karl D. Hammond <karlh@ugcs.caltech.edu>
                          University of Tennessee, Knoxville (USA), 2012
------------------------------------------------------------------------- */

/* This is set of "wrapper" functions to assist LAMMPS.F90, which itself
   provides a (I hope) robust Fortran interface to library.cpp and
   library.h.  All functions herein COULD be added to library.cpp instead of
   including this as a separate file. See the README for instructions. */

#include <mpi.h>
#include "LAMMPS-wrapper2.h"
#include <library.h>
#include <lammps.h>
#include <atom.h>
#include <input.h>
#include <modify.h>
#include <fix.h>
#include <fix_external.h>
#include <compute.h>
#include <modify.h>
#include <error.h>
#include <cstdlib>

using namespace LAMMPS_NS;

extern "C" void f_callback(void *, bigint, int, tagint *, double **, double **);

void lammps_set_callback (void *ptr) {
  class LAMMPS *lmp = (class LAMMPS *) ptr;
  int ifix = lmp->modify->find_fix_by_style("external");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_callback(f_callback, ptr);
  return;
}

void lammps_set_external_vector_length (void *ptr, int n) {
  class LAMMPS *lmp = (class LAMMPS *) ptr;
  int ifix = lmp->modify->find_fix_by_style("external");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_vector_length(n);
  return;
}

void lammps_set_external_vector (void *ptr, int n, double val) {
  class LAMMPS *lmp = (class LAMMPS *) ptr;
  int ifix = lmp->modify->find_fix_by_style("external");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_vector (n, val);
  return;
}

void lammps_set_user_energy (void *ptr, double energy) {
  class LAMMPS *lmp = (class LAMMPS *) ptr;
  int ifix = lmp->modify->find_fix_by_style("external");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_energy_global(energy);
  return;
}

void lammps_set_user_virial (void *ptr, double *virial) {
  class LAMMPS *lmp = (class LAMMPS *) ptr;
  int ifix = lmp->modify->find_fix_by_style("external");
  FixExternal *fix = (FixExternal *) lmp->modify->fix[ifix];
  fix->set_virial_global(virial);
  return;
}

