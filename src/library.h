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

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

#include "mpi.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

void lammps_open(int, char **, MPI_Comm, void **);
void lammps_open_no_mpi(int, char **, void **);
void lammps_close(void *);
void lammps_file(void *, char *);
char *lammps_command(void *, char *);
void lammps_free(void *);

void *lammps_extract_global(void *, char *);
void *lammps_extract_atom(void *, char *);
void *lammps_extract_compute(void *, char *, int, int);
void *lammps_extract_fix(void *, char *, int, int, int, int);
void *lammps_extract_variable(void *, char *, char *);

int lammps_get_natoms(void *);
void lammps_gather_atoms(void *, char *, int, int, void *);
void lammps_scatter_atoms(void *, char *, int, int, void *);

#ifdef __cplusplus
}
#endif

/* ERROR/WARNING messages:

W: Library error in lammps_gather_atoms

This library function cannot be used if atom IDs are not defined
or are not consecutively numbered.

W: Library error in lammps_scatter_atoms

This library function cannot be used if atom IDs are not defined or
are not consecutively numbered, or if no atom map is defined.  See the
atom_modify command for details about atom maps.

*/
