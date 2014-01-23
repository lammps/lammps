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

#include "lmptype.h"
#include "mpi.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "compute.h"
#include "atom.h"
#include "domain.h"
#include "comm.h"
#include "group.h"
#include "atom_masks.h"
#include "memory.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

#define DELTA 4
#define BIG MAXTAGINT

/* ---------------------------------------------------------------------- */

Compute::Compute(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  if (narg < 3) error->all(FLERR,"Illegal compute command");

  // compute ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  for (int i = 0; i < n-1; i++)
    if (!isalnum(id[i]) && id[i] != '_')
      error->all(FLERR,
                 "Compute ID must be alphanumeric or underscore characters");

  igroup = group->find(arg[1]);
  if (igroup == -1) error->all(FLERR,"Could not find compute group ID");
  groupbit = group->bitmask[igroup];

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  // set child class defaults

  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = 0;

  tempflag = pressflag = peflag = 0;
  pressatomflag = peatomflag = 0;
  tempbias = 0;

  timeflag = 0;
  comm_forward = comm_reverse = 0;
  cudable = 0;

  invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_peratom = invoked_local = -1;

  // set modify defaults

  extra_dof = domain->dimension;
  dynamic = 0;

  // setup list of timesteps

  ntime = maxtime = 0;
  tlist = NULL;

  // setup map for molecule IDs

  molmap = NULL;

  datamask = ALL_MASK;
  datamask_ext = ALL_MASK;
}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
  delete [] id;
  delete [] style;

  memory->destroy(tlist);
  memory->destroy(molmap);
}

/* ---------------------------------------------------------------------- */

void Compute::modify_params(int narg, char **arg)
{
  if (narg == 0) error->all(FLERR,"Illegal compute_modify command");

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"extra") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute_modify command");
      extra_dof = force->inumeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"dynamic") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) dynamic = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) dynamic = 1;
      else error->all(FLERR,"Illegal compute_modify command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"thermo") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal compute_modify command");
      if (strcmp(arg[iarg+1],"no") == 0) thermoflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) thermoflag = 1;
      else error->all(FLERR,"Illegal compute_modify command");
      iarg += 2;
    } else error->all(FLERR,"Illegal compute_modify command");
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

void Compute::reset_extra_compute_fix(const char *)
{
  error->all(FLERR,
             "Compute does not allow an extra compute or fix to be reset");
}

/* ----------------------------------------------------------------------
   add ntimestep to list of timesteps the compute will be called on
   do not add if already in list
   search from top downward, since list of times is in decreasing order
------------------------------------------------------------------------- */

void Compute::addstep(bigint ntimestep)
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
    memory->grow(tlist,maxtime,"compute:tlist");
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

int Compute::matchstep(bigint ntimestep)
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

/* ----------------------------------------------------------------------
   identify molecule IDs with atoms in group
   warn if any atom in group has molecule ID = 0
   warn if any molecule has only some atoms in group
   return Ncount = # of molecules with atoms in group
   set molmap to NULL if molecule IDs include all in range from 1 to Ncount
   else: molecule IDs range from idlo to idhi
         set molmap to vector of length idhi-idlo+1
         molmap[id-idlo] = index from 0 to Ncount-1
         return idlo and idhi
------------------------------------------------------------------------- */

int Compute::molecules_in_group(tagint &idlo, tagint &idhi)
{
  int i;

  memory->destroy(molmap);
  molmap = NULL;

  // find lo/hi molecule ID for any atom in group
  // warn if atom in group has ID = 0

  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  tagint lo = BIG;
  tagint hi = -BIG;
  int flag = 0;
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (molecule[i] == 0) flag = 1;
      lo = MIN(lo,molecule[i]);
      hi = MAX(hi,molecule[i]);
    }

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,"Atom with molecule ID = 0 included in "
                   "compute molecule group");

  MPI_Allreduce(&lo,&idlo,1,MPI_LMP_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_LMP_TAGINT,MPI_MAX,world);
  if (idlo == BIG) return 0;

  // molmap = vector of length nlen
  // set to 1 for IDs that appear in group across all procs, else 0

  tagint nlen_tag = idhi-idlo+1;
  if (nlen_tag > MAXSMALLINT) 
    error->all(FLERR,"Too many molecules for compute");
  int nlen = (int) nlen_tag;

  memory->create(molmap,nlen,"compute:molmap");
  for (i = 0; i < nlen; i++) molmap[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      molmap[molecule[i]-idlo] = 1;

  int *molmapall;
  memory->create(molmapall,nlen,"compute:molmapall");
  MPI_Allreduce(molmap,molmapall,nlen,MPI_INT,MPI_MAX,world);

  // nmolecules = # of non-zero IDs in molmap
  // molmap[i] = index of molecule, skipping molecules not in group with -1

  int nmolecules = 0;
  for (i = 0; i < nlen; i++)
    if (molmapall[i]) molmap[i] = nmolecules++;
    else molmap[i] = -1;
  memory->destroy(molmapall);

  // warn if any molecule has some atoms in group and some not in group

  flag = 0;
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) continue;
    if (molecule[i] < idlo || molecule[i] > idhi) continue;
    if (molmap[molecule[i]-idlo] >= 0) flag = 1;
  }

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall && comm->me == 0)
    error->warning(FLERR,
                   "One or more compute molecules has atoms not in group");

  // if molmap simply stores 1 to Nmolecules, then free it

  if (idlo == 1 && idhi == nmolecules && nlen == nmolecules) {
    memory->destroy(molmap);
    molmap = NULL;
  }
  return nmolecules;
}
