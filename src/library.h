/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

// library interface to LAMMPS
// new application-specific functions can be added

#include "mpi.h"

extern "C" {

void lammps_open(int, char **, MPI_Comm); // start LAMMPS w/ command-line args
void lammps_close();                      // shut-down LAMMPS
void lammps_file(char *);                 // execute an input script
char *lammps_command(char *);             // execute a single LAMMPS command

int lammps_get_natoms();               // return # of atoms
void lammps_get_coords(double *);      // retrieve atom coords from all procs
void lammps_put_coords(double *);      // overwrite atom coords on all procs

}
