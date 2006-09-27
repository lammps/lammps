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
#include "minimize.h"
#include "system.h"
#include "domain.h"
#include "update.h"
#include "min.h"
#include "finish.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

void Minimize::command(int narg, char **arg)
{
  if (narg != 3) error->all("Illegal minimize command");

  if (domain->box_exist == 0)
    error->all("Minimize command before simulation box is defined");

  update->tolerance = atof(arg[0]);
  update->nsteps = atoi(arg[1]);
  update->max_eval = atoi(arg[2]);

  update->beginstep = update->firststep = update->ntimestep;
  update->endstep = update->laststep = update->firststep + update->nsteps;

  update->whichflag = 1;

  sys->init();
  update->minimize->run();

  Finish finish;
  finish.end(1);
  update->whichflag = -1;
}
