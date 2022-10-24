/* -----------------------------------------------------------------------
    LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
    https://www.lammps.org/, Sandia National Laboratories
    LAMMPS development team: developers@lammps.org
 
    Copyright (2003) Sandia Corporation.  Under the terms of Contract
    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
    certain rights in this software.  This software is distributed under 
    the GNU General Public License.
 
    See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ------------------------------------------------------------------------
    Contributing author:  Karl D. Hammond <hammondkd@missouri.edu>
                          University of Missouri (USA), 2012
------------------------------------------------------------------------- */

/* This is set of "wrapper" functions to assist LAMMPS.F90, which itself
   provides a (I hope) robust Fortran interface to library.cpp and
   library.h.  All functions herein COULD be added to library.cpp instead of
   including this as a separate file. See the README for instructions. */

#include <mpi.h>
#include "LAMMPS-wrapper.h"
#define LAMMPS_LIB_MPI 1
#include <library.h>
#include <lammps.h>
#include <atom.h>
#include <fix.h>
#include <compute.h>
#include <modify.h>
#include <error.h>
#include <cstdlib>

using namespace LAMMPS_NS;

void lammps_open_fortran_wrapper (int argc, char **argv,
      MPI_Fint communicator, void **ptr)
{
   MPI_Comm C_communicator = MPI_Comm_f2c (communicator);
   lammps_open (argc, argv, C_communicator, ptr);
}

int lammps_get_ntypes (void *ptr)
{
  class LAMMPS *lmp = (class LAMMPS *) ptr;
  int ntypes = lmp->atom->ntypes;
  return ntypes;
}

void lammps_error_all (void *ptr, const char *file, int line, const char *str)
{
   class LAMMPS *lmp = (class LAMMPS *) ptr;
   lmp->error->all (file, line, str);
}

int lammps_extract_compute_vectorsize (void *ptr, char *id, int style)
{
   int *size;
   size = (int *) lammps_extract_compute(ptr, id, style, LMP_SIZE_VECTOR);
   if (size) return *size;
   return 0;
}

void lammps_extract_compute_arraysize (void *ptr, char *id, int style,
      int *nrows, int *ncols)
{
   int *tmp;
   tmp = (int *) lammps_extract_compute(ptr, id, style, LMP_SIZE_ROWS);
   if (tmp) *nrows = *tmp;
   tmp = (int *) lammps_extract_compute(ptr, id, style, LMP_SIZE_COLS);
   if (tmp) *ncols = *tmp;
   return;
}

int lammps_extract_fix_vectorsize (void *ptr, char *id, int style)
{
   int *size;
   size = (int *) lammps_extract_fix(ptr, id, style, LMP_SIZE_VECTOR, 0, 0);
   if (size) return *size;
   return 0;
}

void lammps_extract_fix_arraysize (void *ptr, char *id, int style,
      int *nrows, int *ncols)
{
   int *tmp;
   tmp = (int *) lammps_extract_fix(ptr, id, style, LMP_SIZE_ROWS, 0, 0);
   if (tmp) *nrows = *tmp;
   tmp = (int *) lammps_extract_fix(ptr, id, style, LMP_SIZE_COLS, 0, 0);
   if (tmp) *ncols = *tmp;
   return;
}

/* vim: set ts=3 sts=3 expandtab: */
