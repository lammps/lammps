// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_group.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGroup::FixGroup(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg),
idregion(nullptr), idvar(nullptr), idprop(nullptr)
{
  // dgroupbit = bitmask of dynamic group
  // group ID is last part of fix ID

  auto dgroupid = std::string(id).substr(strlen("GROUP_"));
  gbit = group->bitmask[group->find(dgroupid)];
  gbitinverse = group->inversemask[group->find(dgroupid)];

  // process optional args

  regionflag = 0;
  varflag = 0;
  propflag = 0;
  nevery = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal group command");
      if (domain->find_region(arg[iarg+1]) < 0)
        error->all(FLERR,"Region ID for group dynamic does not exist");
      regionflag = 1;
      delete [] idregion;
      idregion = utils::strdup(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"var") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal group command");
      if (input->variable->find(arg[iarg+1]) < 0)
        error->all(FLERR,"Variable name for group dynamic does not exist");
      varflag = 1;
      delete [] idvar;
      idvar = utils::strdup(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"property") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal group command");
      int flag,cols;
      iprop = atom->find_custom(arg[iarg+1],flag,cols);
      if (iprop < 0 || cols)
        error->all(FLERR,"Custom per-atom vector for group dynamic "
                   "does not exist");
      propflag = 1;
      delete [] idprop;
      idprop = utils::strdup(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal group command");
      nevery = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
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
  delete [] idprop;
}

/* ---------------------------------------------------------------------- */

int FixGroup::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  mask |= POST_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGroup::init()
{
  // parent group cannot be dynamic
  // else order of FixGroup fixes would matter

  if (group->dynamic[igroup])
    error->all(FLERR,"Group dynamic parent group cannot be dynamic");

  if (utils::strmatch(update->integrate_style,"^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  // set current indices for region and variable and custom property

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

  if (propflag) {
    int cols;
    iprop = atom->find_custom(idprop,proptype,cols);
    if (iprop < 0 || cols)
      error->all(FLERR,"Group dynamic command custom property vector does not exist");
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

void FixGroup::setup(int /*vflag*/)
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

void FixGroup::post_integrate_respa(int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa-1) post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixGroup::set_group()
{
  int nlocal = atom->nlocal;

  // invoke atom-style variable if defined
  // set post_integrate flag to 1, then unset after
  // this is for any compute to check if it needs to
  //   operate differently due to invocation this early in timestep
  // e.g. perform ghost comm update due to atoms having just moved

  double *var = nullptr;
  int *ivector = nullptr;
  double *dvector = nullptr;

  if (varflag) {
    update->post_integrate = 1;
    modify->clearstep_compute();
    memory->create(var,nlocal,"fix/group:varvalue");
    input->variable->compute_atom(ivar,igroup,var,1,0);
    modify->addstep_compute(update->ntimestep + nevery);
    update->post_integrate = 0;
  }

  // set ptr to custom atom vector

  if (propflag && !proptype) ivector = atom->ivector[iprop];
  if (propflag && proptype) dvector = atom->dvector[iprop];

  // update region in case it has a variable dependence or is dynamic

  if (regionflag) region->prematch();

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
      if (propflag) {
        if (!proptype && ivector[i] == 0) inflag = 0;
        if (proptype && dvector[i] == 0.0) inflag = 0;
      }
    } else inflag = 0;

    if (inflag) mask[i] |= gbit;
    else mask[i] &= gbitinverse;
  }

  if (varflag) memory->destroy(var);
}

/* ---------------------------------------------------------------------- */

void *FixGroup::extract(const char *str, int &/*unused*/)
{
  if (strcmp(str,"property") == 0 && propflag) return (void *) idprop;
  if (strcmp(str,"variable") == 0 && varflag) return (void *) idvar;
  if (strcmp(str,"region") == 0 && regionflag) return (void *) idregion;
  return nullptr;
}
