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
#include "min.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

Min::Min(LAMMPS *lmp) : Pointers(lmp)
{
  dmax = 0.1;

  elist_atom = NULL;
  vlist_global = vlist_atom = NULL;
}

/* ---------------------------------------------------------------------- */

Min::~Min()
{
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
}

/* ---------------------------------------------------------------------- */

void Min::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal min_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"dmax") == 0) {
      if (iarg+2 > narg) error->all("Illegal min_modify command");
      dmax = atof(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal min_modify command");
  }
}

/* ----------------------------------------------------------------------
   setup lists of computes for global and per-atom PE and pressure
------------------------------------------------------------------------- */

void Min::ev_setup()
{
  delete [] elist_atom;
  delete [] vlist_global;
  delete [] vlist_atom;
  elist_atom = NULL;
  vlist_global = vlist_atom = NULL;

  nelist_atom = 0;
  nvlist_global = nvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peatomflag) nelist_atom++;
    if (modify->compute[i]->pressflag) nvlist_global++;
    if (modify->compute[i]->pressatomflag) nvlist_atom++;
  }

  if (nelist_atom) elist_atom = new Compute*[nelist_atom];
  if (nvlist_global) vlist_global = new Compute*[nvlist_global];
  if (nvlist_atom) vlist_atom = new Compute*[nvlist_atom];

  nelist_atom = 0;
  nvlist_global = nvlist_atom = 0;
  for (int i = 0; i < modify->ncompute; i++) {
    if (modify->compute[i]->peatomflag)
      elist_atom[nelist_atom++] = modify->compute[i];
    if (modify->compute[i]->pressflag)
      vlist_global[nvlist_global++] = modify->compute[i];
    if (modify->compute[i]->pressatomflag)
      vlist_atom[nvlist_atom++] = modify->compute[i];
  }
}

/* ----------------------------------------------------------------------
   set eflag,vflag for current iteration with ntimestep
   always set eflag_global = 1, since need energy every iteration
   eflag = 0 = no energy computation
   eflag = 1 = global energy only
   eflag = 2 = per-atom energy only
   eflag = 3 = both global and per-atom energy
   vflag = 0 = no virial computation (pressure)
   vflag = 1 = global virial with pair portion via sum of pairwise interactions
   vflag = 2 = global virial with pair portion via F dot r including ghosts
   vflag = 4 = per-atom virial only
   vflag = 5 or 6 = both global and per-atom virial
------------------------------------------------------------------------- */

void Min::ev_set(int ntimestep)
{
  int i;

  int eflag_global = 1;

  int eflag_atom = 0;
  for (i = 0; i < nelist_atom; i++)
    if (elist_atom[i]->matchstep(ntimestep)) break;
  if (i < nelist_atom) eflag_atom = 2;

  if (eflag_global) update->eflag_global = update->ntimestep;
  if (eflag_atom) update->eflag_atom = update->ntimestep;
  eflag = eflag_global + eflag_atom;

  int vflag_global = 0;
  for (i = 0; i < nvlist_global; i++)
    if (vlist_global[i]->matchstep(ntimestep)) break;
  if (i < nvlist_global) vflag_global = virial_style;

  int vflag_atom = 0;
  for (i = 0; i < nvlist_atom; i++)
    if (vlist_atom[i]->matchstep(ntimestep)) break;
  if (i < nvlist_atom) vflag_atom = 4;

  if (vflag_global) update->vflag_global = update->ntimestep;
  if (vflag_atom) update->vflag_atom = update->ntimestep;
  vflag = vflag_global + vflag_atom;
}
