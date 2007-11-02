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

#include "stdlib.h"
#include "string.h"
#include "fix_print.h"
#include "update.h"
#include "input.h"
#include "modify.h"
#include "variable.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixPrint::FixPrint(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all("Illegal fix print command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix print command");

  MPI_Comm_rank(world,&me);

  int n = strlen(arg[4]) + 1;
  line = new char[n];
  strcpy(line,arg[4]);

  copy = new char[MAXLINE];
  work = new char[MAXLINE];
}

/* ---------------------------------------------------------------------- */

FixPrint::~FixPrint()
{
  delete [] line;
  delete [] copy;
  delete [] work;
}

/* ---------------------------------------------------------------------- */

int FixPrint::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPrint::end_of_step()
{
  // make a copy of line to work on
  // substitute for $ variables (no printing)
  // append a newline and print final copy
  // variable evaluation may invoke a compute that affects Verlet::eflag,vflag

  modify->clearstep_compute();

  strcpy(copy,line);
  input->substitute(copy,0);
  strcat(copy,"\n");

  modify->addstep_compute(update->ntimestep + nevery);

  if (me == 0) {
    if (screen) fprintf(screen,copy);
    if (logfile) fprintf(logfile,copy);
  }
}
