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
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "region.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGroup::FixGroup(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idregion(nullptr), idvar(nullptr), idprop(nullptr), idexclude(nullptr),
    region(nullptr), list(nullptr)
{
  // gbit = bitmask of dynamic group
  // group ID is last part of fix ID

  auto dgroupid = std::string(id).substr(strlen("GROUP_"));
  gbit = group->get_bitmask_by_id(FLERR, dgroupid, "dynamic group");
  gbitinverse = group->get_inversemask_by_id(FLERR, dgroupid, "dynamic group");

  comm_forward = 1;

  // process optional args

  regionflag = 0;
  varflag = 0;
  propflag = 0;
  moleculeflag = 0;
  withinflag = 0;
  excludeflag = 0;
  cutoff = 0.0;
  nevery = 1;
  excludebit = 0;

  int ioffset = 0;
  if (lmp->input->arg) {
    for (int i = 0; i < lmp->input->narg; ++i)
      if (lmp->input->arg[i] == arg[3]) ioffset = i;
  }

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic region", error);
      if (!domain->get_region_by_id(arg[iarg + 1]))
        error->all(FLERR, ioffset + iarg + 1, "Region {} for dynamic group {} does not exist",
                   arg[iarg + 1], dgroupid);
      regionflag = 1;
      delete[] idregion;
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "var") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic var", error);
      if (input->variable->find(arg[iarg + 1]) < 0)
        error->all(FLERR, ioffset + iarg + 1, "Variable '{}' for dynamic group {} does not exist",
                   arg[iarg + 1], dgroupid);
      varflag = 1;
      delete[] idvar;
      idvar = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "property") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic property", error);
      int flag, cols;
      iprop = atom->find_custom(arg[iarg + 1], flag, cols);
      if (iprop < 0 || cols)
        error->all(FLERR, ioffset + iarg + 1,
                   "Custom per-atom vector {} for dynamic group {} does not exist", arg[iarg + 1],
                   dgroupid);
      propflag = 1;
      delete[] idprop;
      idprop = utils::strdup(arg[iarg + 1]);
      iarg += 2;

    } else if (strcmp(arg[iarg], "every") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic every", error);
      nevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (nevery <= 0)
        error->all(FLERR, ioffset + iarg + 1, "Illegal every value {} for dynamic group {}", nevery,
                   dgroupid);
      iarg += 2;
    } else if (strcmp(arg[iarg], "include") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic include", error);
      if (strcmp(arg[iarg + 1], "molecule") == 0) {
        moleculeflag = 1;
        if (!atom->molecule_flag)
          error->all(FLERR, ioffset + iarg + 1,
                     "Dynamic Group include molecule setting requires atom attribute molecule");
        iarg += 2;
      } else {
        error->all(FLERR, ioffset + iarg + 1, "Unknown include setting {} in dynamic group command",
                   arg[iarg + 1]);
      }
    } else if (strcmp(arg[iarg], "within") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic within", error);
      withinflag = 1;
      cutoff = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (cutoff <= 0.0)
        error->all(FLERR, ioffset + iarg + 1, "Illegal within cutoff value {} for dynamic group {}",
                   cutoff, dgroupid);
      iarg += 2;
    } else if (strcmp(arg[iarg], "exclude") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "group dynamic exclude", error);
      excludeflag = 1;
      delete[] idexclude;
      idexclude = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, ioffset + iarg, "Unknown keyword {} in dynamic group command", arg[iarg]);
  }
}

/* ---------------------------------------------------------------------- */

FixGroup::~FixGroup()
{
  delete[] idregion;
  delete[] idvar;
  delete[] idprop;
  delete[] idexclude;
}

/* ---------------------------------------------------------------------- */

int FixGroup::setmask()
{
  int mask = 0;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGroup::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
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

  if (withinflag) {
    NeighRequest *req = nullptr;
    if (nevery == 1) {
      req = neighbor->add_request(this, NeighConst::REQ_FULL);
    } else {
      req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
    }
    req->set_cutoff(cutoff);
  }

  if (excludeflag) excludebit = group->get_bitmask_by_id(FLERR, idexclude, "group dynamic exclude");
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

  // re-build occasional neighbor list for "within" processing

  if (withinflag && (nevery != 1)) neighbor->build_one(list);

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

  // ensure ghost atom masks are also updated

  comm->forward_comm(this);

  // select additional atoms that are within cutoff distance of already selected atoms

  if (withinflag) {
    int i, j, ii, jj, inum, jnum;
    const int *ilist, *jlist;
    double dx, dy, dz, rsq;

    inum = list->inum;
    ilist = list->ilist;
    const int *const numneigh = list->numneigh;
    const int *const *const firstneigh = list->firstneigh;

    const double cutsq = cutoff * cutoff;
    for (ii = 0; ii < inum; ++ii) {
      i = ilist[ii];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      for (jj = 0; jj < jnum; ++jj) {
        j = jlist[jj] & NEIGHMASK;

        // check pairs where only one of the two atoms is in the parent group and passes the
        // selection criteria from above.  we cannot use mask[] & gbit for this for local atoms,
        // since we add atoms to the group and we don't want to add neighbors of these added atoms
        if ((mask[i] & groupbit) && !(mask[j] & groupbit)) {
          inflag = 1;
          if (regionflag && !region->match(x[i][0], x[i][1], x[i][2])) inflag = 0;
          if (varflag && var[i] == 0.0) inflag = 0;
          if (propflag) {
            if (!proptype && ivector[i] == 0) inflag = 0;
            if (proptype && dvector[i] == 0.0) inflag = 0;
          }
          if (inflag && (j < nlocal)) {
            dx = x[j][0] - x[i][0];
            dy = x[j][1] - x[i][1];
            dz = x[j][2] - x[i][2];
            rsq = dx * dx + dy * dy + dz * dz;
            if (rsq <= cutsq) mask[j] |= gbit;
          }
        } else if (!(mask[i] & groupbit) && (mask[j] & groupbit)) {
          inflag = 1;
          // for ghost atoms, we do not have access to all data, but we can use gbit instead
          // since it has already been forward communicated and we will only select local atoms
          if (j >= nlocal) {
            if (!(mask[j] & gbit)) inflag = 0;
          } else {
            if (regionflag && !region->match(x[j][0], x[j][1], x[j][2])) inflag = 0;
            if (varflag && var[j] == 0.0) inflag = 0;
            if (propflag) {
              if (!proptype && ivector[j] == 0) inflag = 0;
              if (proptype && dvector[j] == 0.0) inflag = 0;
            }
          }
          if (inflag) {
            dx = x[j][0] - x[i][0];
            dy = x[j][1] - x[i][1];
            dz = x[j][2] - x[i][2];
            rsq = dx * dx + dy * dy + dz * dz;
            if (rsq <= cutsq) mask[i] |= gbit;
          }
        }
      }
    }

    // we need a second forward communication, since we could only update the masks of local atoms
    comm->forward_comm(this);
  }

  // no longer needed

  if (varflag) memory->destroy(var);

  // add atoms that have the same molecule ID as selected atoms
  if (moleculeflag) group->add_molecules(0, gbit);

  // exclude selected atoms in the excluded group
  if (excludeflag) {
    int nall = nlocal + atom->nghost;
    for (int i = 0; i < nall; ++i) {
      if (mask[i] & gbit) {
        if (mask[i] & excludebit) mask[i] &= gbitinverse;
      }
    }
  }
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
