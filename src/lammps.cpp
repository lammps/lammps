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

// parent class for all of LAMMPS
// all other classes inherit from this
// contains static ptrs to single instance of other classes
// contains MPI communicator and file handles for my world of procs

#include "mpi.h"
#include "stdlib.h"
#include "lammps.h"
#include "system.h"

// set static ptrs to NULL

System *LAMMPS::sys = NULL;
Universe *LAMMPS::universe = NULL;
Input *LAMMPS::input = NULL;
Memory *LAMMPS::memory = NULL;
Error *LAMMPS::error = NULL;

Atom *LAMMPS::atom = NULL;
Update *LAMMPS::update = NULL;
Neighbor *LAMMPS::neighbor = NULL;
Comm *LAMMPS::comm = NULL;
Domain *LAMMPS::domain = NULL;
Force *LAMMPS::force = NULL;
Modify *LAMMPS::modify = NULL;
Group *LAMMPS::group = NULL;
Output *LAMMPS::output = NULL;
Timer *LAMMPS::timer = NULL;

MPI_Comm LAMMPS::world = 0;
FILE *LAMMPS::infile = NULL;
FILE *LAMMPS::screen = NULL;
FILE *LAMMPS::logfile = NULL;

/* ---------------------------------------------------------------------- */

void LAMMPS::open(int narg, char **arg, MPI_Comm communicator)
{
  sys = new System();
  sys->open(narg,arg,communicator);
  sys->create();
}

/* ---------------------------------------------------------------------- */

void LAMMPS::close()
{
  sys->destroy();
  sys->close();
  delete sys;
}
