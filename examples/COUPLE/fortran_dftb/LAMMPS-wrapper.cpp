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
#include "LAMMPS-wrapper.h"
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
   class LAMMPS *lmp = (class LAMMPS *) ptr;
   int icompute = lmp->modify->find_compute(id);
   if ( icompute < 0 ) return 0;
   class Compute *compute = lmp->modify->compute[icompute];

   if ( style == 0 )
   {
      if ( !compute->vector_flag )
         return 0;
      else
         return compute->size_vector;
   }
   else if ( style == 1 )
   {
      return lammps_get_natoms (ptr);
   }
   else if ( style == 2 )
   {
      if ( !compute->local_flag )
         return 0;
      else
         return compute->size_local_rows;
   }
   else
      return 0;
}

void lammps_extract_compute_arraysize (void *ptr, char *id, int style,
      int *nrows, int *ncols)
{
   class LAMMPS *lmp = (class LAMMPS *) ptr;
   int icompute = lmp->modify->find_compute(id);
   if ( icompute < 0 )
   {
      *nrows = 0;
      *ncols = 0;
   }
   class Compute *compute = lmp->modify->compute[icompute];

   if ( style == 0 )
   {
      if ( !compute->array_flag )
      {
         *nrows = 0;
         *ncols = 0;
      }
      else
      {
         *nrows = compute->size_array_rows;
         *ncols = compute->size_array_cols;
      }
   }
   else if ( style == 1 )
   {
      if ( !compute->peratom_flag )
      {
         *nrows = 0;
         *ncols = 0;
      }
      else
      {
         *nrows = lammps_get_natoms (ptr);
         *ncols = compute->size_peratom_cols;
      }
   }
   else if ( style == 2 )
   {
      if ( !compute->local_flag )
      {
         *nrows = 0;
         *ncols = 0;
      }
      else
      {
         *nrows = compute->size_local_rows;
         *ncols = compute->size_local_cols;
      }
   }
   else
   {
      *nrows = 0;
      *ncols = 0;
   }

   return;
}

int lammps_extract_fix_vectorsize (void *ptr, char *id, int style)
{
   class LAMMPS *lmp = (class LAMMPS *) ptr;
   int ifix = lmp->modify->find_fix(id);
   if ( ifix < 0 ) return 0;
   class Fix *fix = lmp->modify->fix[ifix];

   if ( style == 0 )
   {
      if ( !fix->vector_flag )
         return 0;
      else
         return fix->size_vector;
   }
   else if ( style == 1 )
   {
      return lammps_get_natoms (ptr);
   }
   else if ( style == 2 )
   {
      if ( !fix->local_flag )
         return 0;
      else
         return fix->size_local_rows;
   }
   else
      return 0;
}

void lammps_extract_fix_arraysize (void *ptr, char *id, int style,
      int *nrows, int *ncols)
{
   class LAMMPS *lmp = (class LAMMPS *) ptr;
   int ifix = lmp->modify->find_fix(id);
   if ( ifix < 0 )
   {
      *nrows = 0;
      *ncols = 0;
   }
   class Fix *fix = lmp->modify->fix[ifix];

   if ( style == 0 )
   {
      if ( !fix->array_flag )
      {
         *nrows = 0;
         *ncols = 0;
      }
      else
      {
         *nrows = fix->size_array_rows;
         *ncols = fix->size_array_cols;
      }
   }
   else if ( style == 1 )
   {
      if ( !fix->peratom_flag )
      {
         *nrows = 0;
         *ncols = 0;
      }
      else
      {
         *nrows = lammps_get_natoms (ptr);
         *ncols = fix->size_peratom_cols;
      }
   }
   else if ( style == 2 )
   {
      if ( !fix->local_flag )
      {
         *nrows = 0;
         *ncols = 0;
      }
      else
      {
         *nrows = fix->size_local_rows;
         *ncols = fix->size_local_cols;
      }
   }
   else
   {
      *nrows = 0;
      *ncols = 0;
   }

   return;

}

/* vim: set ts=3 sts=3 expandtab: */
