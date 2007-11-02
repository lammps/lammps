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
#include "integrate.h"
#include "compute.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Integrate::Integrate(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  elist = vlist = NULL;
}

/* ---------------------------------------------------------------------- */

Integrate::~Integrate()
{
  delete [] elist;
  delete [] vlist;
}

/* ----------------------------------------------------------------------
   set eflag if a pe compute is called this timestep
   set vflag if a pressure compute is called this timestep
------------------------------------------------------------------------- */

void Integrate::ev_set(int ntimestep)
{
  int i;

  eflag = 0;
  for (i = 0; i < nelist; i++)
    if (elist[i]->match_step(ntimestep)) break;
  if (i < nelist) eflag = 1;

  vflag = 0;
  for (i = 0; i < nvlist; i++)
    if (vlist[i]->match_step(ntimestep)) break;
  if (i < nvlist) vflag = virial_style;
}
