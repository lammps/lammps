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

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "temperature.h"
#include "atom.h"
#include "group.h"
#include "modify.h"
#include "domain.h"
#include "lattice.h"
#include "fix.h"
#include "error.h"

/* ---------------------------------------------------------------------- */

Temperature::Temperature(int narg, char **arg)
{
  if (narg < 3) error->all("Illegal temperature command");

  // temperature ID, group, and style

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = group->find(arg[1]);
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  // set modify defaults

  extra_dof = 3;

  // set input line defaults

  scaleflag = 1;

  // read options from end of input line

  if (strcmp(style,"full") == 0) options(narg-3,&arg[3]);
  else if (strcmp(style,"partial") == 0) options(narg-6,&arg[6]);
  else if (strcmp(style,"ramp") == 0) options(narg-9,&arg[9]);
  else if (strcmp(style,"region") == 0) options(narg-4,&arg[4]);

  // set scaling for RAMP style

  if (strcmp(style,"ramp") == 0) {

    if (scaleflag && domain->lattice == NULL)
      error->all("Use of temperature ramp with undefined lattice");

    if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
    }
    else xscale = yscale = zscale = 1.0;
  }
}

/* ---------------------------------------------------------------------- */

Temperature::~Temperature()
{
  delete [] id;
  delete [] style;
}

/* ---------------------------------------------------------------------- */

void Temperature::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal temp_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"extra") == 0) {
      if (iarg+2 > narg) error->all("Illegal temp_modify command");
      extra_dof = atoi(arg[iarg+1]);
      iarg += 2;
    } else error->all("Illegal temp_modify command");
  }
}

/* ----------------------------------------------------------------------
   ncount = atoms in Temperature group 
------------------------------------------------------------------------- */

void Temperature::count_atoms()
{
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int icount = 0;
  for (int i = 0; i < nlocal; i++) if (mask[i] & groupbit) icount++;
  double rcount = icount;
  MPI_Allreduce(&rcount,&ncount,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ----------------------------------------------------------------------
   count degrees of freedom subtracted by fixes
------------------------------------------------------------------------- */

void Temperature::count_fix()
{
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of temperature input line 
------------------------------------------------------------------------- */

void Temperature::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal temperature command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal temperature command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal temperature command");
      iarg += 2;
    } else error->all("Illegal temperature command");
  }
}
