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

#include "stdlib.h"
#include "string.h"
#include "fix_print.h"
#include "update.h"
#include "input.h"
#include "variable.h"
#include "error.h"

#define MAXLINE 1024

/* ---------------------------------------------------------------------- */

FixPrint::FixPrint(int narg, char **arg) : Fix(narg, arg)
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

  strcpy(copy,line);
  input->substitute(copy,0);
  strcat(copy,"\n");

  if (me == 0) {
    if (screen) fprintf(screen,copy);
    if (logfile) fprintf(logfile,copy);
  }
}
