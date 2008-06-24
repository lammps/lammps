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

/* 
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

#include "mpi.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

void lammps_open(int, char **, MPI_Comm, void **);  /* start-up LAMMPS */
void lammps_close(void *);                          /* shut-down LAMMPS */
void lammps_file(void *, char *);                   /* run an input script */
char *lammps_command(void *, char *);               /* execute a command */

int lammps_get_natoms(void *);              /* return # of atoms */
void lammps_get_coords(void *, double *);   /* get atom x from all procs */
void lammps_put_coords(void *, double *);   /* put atom x on all procs */

#ifdef __cplusplus
}
#endif

