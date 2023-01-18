/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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

FixGroup::FixGroup(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idregion(nullptr), idvar(nullptr), idprop(nullptr), region(nullptr)
{
  // dgroupbit = bitmask of dynamic group
  // group ID is last part of fix ID

  auto dgroupid = std::string(id).substr(strlen("GROUP_"));
  gbit = group->bitmask[group->find(dgroupid)];
  gbitinverse = group->inversemask[group->find(dgroupid)];

  comm_forward = 1;

  // process optional args

  regionflag = 0;
  varflag = 0;
  propflag = 0;
  nevery = 1;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic region", error);
      if (!domain->get_region_by_id(arg[iarg + 1]))
        error->all(FLERR, "Region {} for dynamic group {} does not exist", arg[iarg + 1], dgroupid);
      regionflag = 1;
      delete[] idregion;
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "var") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic var", error);
      if (input->variable->find(arg[iarg + 1]) < 0)
        error->all(FLERR, "Variable '{}' for dynamic group {} does not exist", arg[iarg + 1],
                   dgroupid);
      varflag = 1;
      delete[] idvar;
      idvar = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "property") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic property", error);
      int flag, cols;
      iprop = atom->find_custom(arg[iarg + 1], flag, cols);
      if (iprop < 0 || cols)
        error->all(FLERR, "Custom per-atom vector {} for dynamic group {} does not exist",
                   arg[iarg + 1], dgroupid);
      propflag = 1;
      delete[] idprop;
      idprop = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "every") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic every", error);
      nevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (nevery <= 0)
        error->all(FLERR, "Illegal every value {} for dynamic group {}", nevery, dgroupid);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown keyword {} in dynamic group command", arg[iarg]);
  }
}

/* ---------------------------------------------------------------------- */

FixGroup::~FixGroup()
{
  delete[] idregion;
  delete[] idvar;
  delete[] idprop;
}

/* ---------------------------------------------------------------------- */

int FixGroup::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGroup::init()
{
  std::string dyngroup = group->names[igroup];
  // parent group cannot be dynamic
  // else order of FixGroup fixes would matter

  if (group->dynamic[igroup])
    error->all(FLERR, "Dynamic group parent group {} cannot be dynamic", dyngroup);

  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;

  // set current indices for region and variable and custom property

  if (regionflag) {
    region = domain->get_region_by_id(idregion);
    if (!region)
      error->all(FLERR, "Region {} for dynamic group {} does not exist", idregion, dyngroup);
  }

  if (varflag) {
    ivar = input->variable->find(idvar);
    if (ivar < 0)
      error->all(FLERR, "Variable '{}' for dynamic group {} does not exist", idvar, dyngroup);
    if (!input->variable->atomstyle(ivar))
      error->all(FLERR, "Variable '{}' for dynamic group is of incompatible style");
  }

  if (propflag) {
    int cols;
    iprop = atom->find_custom(idprop, proptype, cols);
    if (iprop < 0 || cols)
      error->all(FLERR, "Custom per-atom property vector {} for dynamic group {} does not exist",
                 idprop, dyngroup);
  }
}

/* ----------------------------------------------------------------------
   assign atoms to group
------------------------------------------------------------------------- */

void FixGroup::setup(int /*vflag*/)
{
  set_group();
}

/* ---------------------------------------------------------------------- */

void FixGroup::post_force(int /*vflag*/)
{
  // only assign atoms to group on steps that are multiples of nevery

  if (update->ntimestep % nevery == 0) set_group();
}

/* ---------------------------------------------------------------------- */

void FixGroup::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGroup::set_group()
{
  int nlocal = atom->nlocal;

  // invoke atom-style variable if defined
  // NOTE: after variable invocation could reset invoked computes to not-invoked
  //   this would avoid an issue where other post-force fixes
  //   change the compute result since it will not be re-invoked at end-of-step,
  //   e.g. if compute pe/atom includes pe contributions from fixes

  double *var = nullptr;
  int *ivector = nullptr;
  double *dvector = nullptr;

  if (varflag) {
    modify->clearstep_compute();
    memory->create(var, nlocal, "fix/group:var");
    input->variable->compute_atom(ivar, igroup, var, 1, 0);
    modify->addstep_compute(update->ntimestep + nevery);
  }

  // set ptr to custom atom vector

  if (propflag && !proptype) ivector = atom->ivector[iprop];
  if (propflag && proptype) dvector = atom->dvector[iprop];

  // update region in case it has a variable dependence or is dynamic

  if (regionflag) region->prematch();

  // set mask for each atom
  // only in group if in parent group, in region, variable is non-zero

  double **x = atom->x;
  int *mask = atom->mask;
  int inflag;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      inflag = 1;
      if (regionflag && !region->match(x[i][0], x[i][1], x[i][2])) inflag = 0;
      if (varflag && var[i] == 0.0) inflag = 0;
      if (propflag) {
        if (!proptype && ivector[i] == 0) inflag = 0;
        if (proptype && dvector[i] == 0.0) inflag = 0;
      }
    } else
      inflag = 0;

    if (inflag)
      mask[i] |= gbit;
    else
      mask[i] &= gbitinverse;
  }

  if (varflag) memory->destroy(var);

  // insure ghost atom masks are also updated

  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int FixGroup::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  int *mask = atom->mask;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(mask[j]).d;
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixGroup::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;

  m = 0;
  last = first + n;

  int *mask = atom->mask;

  for (i = first; i < last; i++) mask[i] = (int) ubuf(buf[m++]).i;
}

/* ---------------------------------------------------------------------- */

void *FixGroup::extract(const char *str, int & /*unused*/)
{
  if (strcmp(str, "property") == 0 && propflag) return (void *) idprop;
  if (strcmp(str, "variable") == 0 && varflag) return (void *) idvar;
  if (strcmp(str, "region") == 0 && regionflag) return (void *) idregion;
  return nullptr;
}
