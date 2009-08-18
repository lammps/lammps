/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory. // 
------------------------------------------------------------------------- */

#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "lattice.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

Compute::Compute(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  if (narg < 3) error->all("Illegal compute command");

  // compute ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all("Compute ID must be alphanumeric or underscore characters");

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all("Could not find compute group ID");
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  // set child class defaults

  vector = NULL;
  scalar_atom = NULL;
  vector_atom = NULL;
  
  scalar_flag = vector_flag = peratom_flag = 0;
  tempflag = pressflag = peflag = 0;
  pressatomflag = peatomflag = 0;
  tempbias = 0;

  timeflag = 0;
  comm_forward = comm_reverse = 0;

  invoked_scalar = invoked_vector = invoked_peratom = -1;

  // set modify defaults

  extra_dof = domain->dimension;
  dynamic = 0;

  // setup list of timesteps
  
  ntime = maxtime = 0;
  tlist = NULL;
}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
  delete [] id;
  delete [] style;

  memory->sfree(tlist);
}

/* ---------------------------------------------------------------------- */

void Compute::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all("Illegal compute_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"extra") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute_modify command");
      extra_dof = atoi(arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dynamic") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) dynamic = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) dynamic = 1;
      else error->all("Illegal compute_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"thermo") == 0) {
      if (iarg+2 > narg) error->all("Illegal compute_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) thermoflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) thermoflag = 1;
      else error->all("Illegal compute_modify command");
      iarg += 2;
    } else error->all("Illegal compute_modify command");
  }
}

/* ----------------------------------------------------------------------
   reset extra_dof to its default value
------------------------------------------------------------------------- */

void Compute::reset_extra_dof()
{
  extra_dof = domain->dimension;
}

/* ---------------------------------------------------------------------- */

void Compute::reset_extra_compute_fix(char *)
{
  error->all("Compute does not allow an extra compute or fix to be reset");
}

/* ----------------------------------------------------------------------
   add ntimestep to list of timesteps the compute will be called on
   do not add if already in list
   search from top downward, since list of times is in decreasing order
------------------------------------------------------------------------- */

void Compute::addstep(int ntimestep)
{
  // i = location in list to insert ntimestep

  int i;
  for (i = ntime-1; i >= 0; i--) {
    if (ntimestep == tlist[i]) return;
    if (ntimestep < tlist[i]) break;
  }
  i++;

  // extend list as needed

  if (ntime == maxtime) {
    maxtime += DELTA;
    tlist = (int *)
      memory->srealloc(tlist,maxtime*sizeof(int),"compute:tlist");
  }

  // move remainder of list upward and insert ntimestep

  for (int j = ntime-1; j >= i; j--) tlist[j+1] = tlist[j];
  tlist[i] = ntimestep;
  ntime++;
}

/* ----------------------------------------------------------------------
   return 1/0 if ntimestep is or is not in list of calling timesteps
   if value(s) on top of list are less than ntimestep, delete them
   search from top downward, since list of times is in decreasing order
------------------------------------------------------------------------- */

int Compute::matchstep(int ntimestep)
{
  for (int i = ntime-1; i >= 0; i--) {
    if (ntimestep < tlist[i]) return 0;
    if (ntimestep == tlist[i]) return 1;
    if (ntimestep > tlist[i]) ntime--;
  }
  return 0;
}

/* ----------------------------------------------------------------------
   clean out list of timesteps to call the compute on
------------------------------------------------------------------------- */

void Compute::clearstep()
{
  ntime = 0;
}
