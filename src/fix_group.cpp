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

#include "string.h"
#include "fix_group.h"
#include "group.h"
#include "update.h"
#include "atom.h"
#include "force.h"
#include "comm.h"
#include "domain.h"
#include "region.h"
#include "modify.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGroup::FixGroup(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  // dgroupbit = bitmask of dynamic group
  // group ID is last part of fix ID

  int n = strlen(id) - strlen("GROUP_") + 1;
  char *dgroup = new char[n];
  strcpy(dgroup,&id[strlen("GROUP_")]);
  gbit = group->bitmask[group->find(dgroup)];
  gbitinverse = group->inversemask[group->find(dgroup)];
  delete [] dgroup;

  // process optional args

  regionflag = 0;
  idregion = NULL;
  varflag = 0;
  idvar = NULL;
  nevery = 1;
  
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal group command");
      if (domain->find_region(arg[iarg+1]) < 0)
        error->all(FLERR,"Region ID for group dynamic does not exist");
      regionflag = 1;
      delete [] idregion;
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"var") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal group command");
      if (input->variable->find(arg[iarg+1]) < 0)
        error->all(FLERR,"Variable name for group dynamic does not exist");
      varflag = 1;
      delete [] idvar;
      int n = strlen(arg[iarg+1]) + 1;
      idvar = new char[n];
      strcpy(idvar,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal group command");
      nevery = force->inumeric(FLERR,arg[iarg+1]);
      if (nevery <= 0) error->all(FLERR,"Illegal group command");
      iarg += 2;
    } else error->all(FLERR,"Illegal group command");
  }
}

/* ---------------------------------------------------------------------- */

FixGroup::~FixGroup()
{
  delete [] idregion;
  delete [] idvar;
}

/* ---------------------------------------------------------------------- */

int FixGroup::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGroup::init()
{
  // parent group cannot be dynamic
  // else order of FixGroup fixes would matter

  if (group->dynamic[igroup])
    error->all(FLERR,"Group dynamic parent group cannot be dynamic");

  // set current indices for region and variable

  if (regionflag) {
    iregion = domain->find_region(idregion);
    if (iregion < 0)
      error->all(FLERR,"Region ID for group dynamic does not exist");
    region = domain->regions[iregion];
  }

  if (varflag) {
    ivar = input->variable->find(idvar);
    if (ivar < 0)
      error->all(FLERR,"Variable name for group dynamic does not exist");
    if (!input->variable->atomstyle(ivar))
      error->all(FLERR,"Variable for group dynamic is invalid style");
  }

  // warn if any FixGroup is not at tail end of all post_integrate fixes

  Fix **fix = modify->fix;
  int *fmask = modify->fmask;
  int nfix = modify->nfix;

  int n = 0;
  for (int i = 0; i < nfix; i++) if (POST_INTEGRATE & fmask[i]) n++;
  int warn = 0;
  for (int i = 0; i < nfix; i++) {
    if (POST_INTEGRATE & fmask[i]) {
      for (int j = i+1; j < nfix; j++) {
        if (POST_INTEGRATE & fmask[j]) {
          if (strstr(fix[j]->id,"GROUP_") != fix[j]->id) warn = 1;
        }
      }
    }
  }

  if (warn && comm->me == 0) 
    error->warning(FLERR,"One or more dynamic groups may not be "
                   "updated at correct point in timestep");
}

/* ----------------------------------------------------------------------
   assign atoms to group
------------------------------------------------------------------------- */

void FixGroup::setup(int vflag)
{
  set_group();
}

/* ---------------------------------------------------------------------- */

void FixGroup::post_integrate()
{
  // only assign atoms to group on steps that are multiples of nevery

  if (update->ntimestep % nevery == 0) set_group();
}

/* ---------------------------------------------------------------------- */

void FixGroup::set_group()
{
  int nlocal = atom->nlocal;

  // invoke atom-style variable if defined

  double *var = NULL;

  if (varflag) {
    modify->clearstep_compute();
    memory->create(var,nlocal,"fix/group:varvalue");
    input->variable->compute_atom(ivar,igroup,var,1,0);
    modify->addstep_compute(update->ntimestep + nevery);
  }
  
  // set mask for each atom
  // only in group if in parent group, in region, variable is non-zero
  // if compute, fix, etc needs updated masks of ghost atoms,
  // it must do forward_comm() to update them

  double **x = atom->x;
  int *mask = atom->mask;
  int inflag;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      inflag = 1;
      if (regionflag && !region->match(x[i][0],x[i][1],x[i][2])) inflag = 0;
      if (varflag && var[i] == 0.0) inflag = 0;
    } else inflag = 0;

    if (inflag) mask[i] |= gbit;
    else mask[i] &= gbitinverse;
  }

  if (varflag) memory->destroy(var);
}
