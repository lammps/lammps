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

#include "fix_setforce.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
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

enum { NONE, CONSTANT, EQUAL, ATOM };

/* ---------------------------------------------------------------------- */

FixSetForce::FixSetForce(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr), idregion(nullptr),
    region(nullptr), sforce(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix setforce", error);

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = nlevels_respa = 0;

  if (utils::strmatch(arg[3], "^v_")) {
    xstr = utils::strdup(arg[3] + 2);
  } else if (strcmp(arg[3], "NULL") == 0) {
    xstyle = NONE;
  } else {
    xvalue = utils::numeric(FLERR, arg[3], false, lmp);
    xstyle = CONSTANT;
  }
  if (utils::strmatch(arg[4], "^v_")) {
    ystr = utils::strdup(arg[4] + 2);
  } else if (strcmp(arg[4], "NULL") == 0) {
    ystyle = NONE;
  } else {
    yvalue = utils::numeric(FLERR, arg[4], false, lmp);
    ystyle = CONSTANT;
  }
  if (utils::strmatch(arg[5], "^v_")) {
    zstr = utils::strdup(arg[5] + 2);
  } else if (strcmp(arg[5], "NULL") == 0) {
    zstyle = NONE;
  } else {
    zvalue = utils::numeric(FLERR, arg[5], false, lmp);
    zstyle = CONSTANT;
  }

  // optional args

  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix setforce region", error);
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix setforce does not exist", arg[iarg + 1]);
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown fix setforce keyword: {}", arg[iarg]);
  }

  force_flag = 0;
  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;

  maxatom = 1;
  memory->create(sforce, maxatom, 3, "setforce:sforce");
}

/* ---------------------------------------------------------------------- */

FixSetForce::~FixSetForce()
{
  if (copymode) return;

  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] idregion;
  memory->destroy(sforce);
}

/* ---------------------------------------------------------------------- */

int FixSetForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixSetForce::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR, "Variable {} for fix setforce does not exist", xstr);
    if (input->variable->equalstyle(xvar))
      xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar))
      xstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix setforce is invalid style", xstr);
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR, "Variable {} for fix setforce does not exist", ystr);
    if (input->variable->equalstyle(yvar))
      ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar))
      ystyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix setforce is invalid style", ystr);
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR, "Variable {} for fix setforce does not exist", zstr);
    if (input->variable->equalstyle(zvar))
      zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar))
      zstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix setforce is invalid style", zstr);
  }

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix setforce does not exist", idregion);
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else
    varflag = CONSTANT;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
    if (respa_level >= 0)
      ilevel_respa = MIN(respa_level, nlevels_respa - 1);
    else
      ilevel_respa = nlevels_respa - 1;
  }

  // cannot use non-zero forces for a minimization since no energy is integrated
  // use fix addforce instead

  int flag = 0;
  if (update->whichflag == 2) {
    if (xstyle == EQUAL || xstyle == ATOM) flag = 1;
    if (ystyle == EQUAL || ystyle == ATOM) flag = 1;
    if (zstyle == EQUAL || zstyle == ATOM) flag = 1;
    if (xstyle == CONSTANT && xvalue != 0.0) flag = 1;
    if (ystyle == CONSTANT && yvalue != 0.0) flag = 1;
    if (zstyle == CONSTANT && zvalue != 0.0) flag = 1;
  }
  if (flag) error->all(FLERR, "Cannot use non-zero forces in an energy minimization");
}

/* ---------------------------------------------------------------------- */

void FixSetForce::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel);
      post_force_respa(vflag, ilevel, 0);
      (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel);
    }
}

/* ---------------------------------------------------------------------- */

void FixSetForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixSetForce::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // update region if necessary

  if (region) region->prematch();

  // reallocate sforce array if necessary

  if (varflag == ATOM && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(sforce);
    memory->create(sforce, maxatom, 3, "setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        foriginal[0] += f[i][0];
        foriginal[1] += f[i][1];
        foriginal[2] += f[i][2];
        if (xstyle) f[i][0] = xvalue;
        if (ystyle) f[i][1] = yvalue;
        if (zstyle) f[i][2] = zvalue;
      }

    // variable force, wrap with clear/add

  } else {

    modify->clearstep_compute();

    if (xstyle == EQUAL)
      xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar, igroup, &sforce[0][0], 3, 0);
    if (ystyle == EQUAL)
      yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar, igroup, &sforce[0][1], 3, 0);
    if (zstyle == EQUAL)
      zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar, igroup, &sforce[0][2], 3, 0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        foriginal[0] += f[i][0];
        foriginal[1] += f[i][1];
        foriginal[2] += f[i][2];
        if (xstyle == ATOM)
          f[i][0] = sforce[i][0];
        else if (xstyle)
          f[i][0] = xvalue;
        if (ystyle == ATOM)
          f[i][1] = sforce[i][1];
        else if (ystyle)
          f[i][1] = yvalue;
        if (zstyle == ATOM)
          f[i][2] = sforce[i][2];
        else if (zstyle)
          f[i][2] = zvalue;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetForce::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  // set force to desired value on requested level, 0.0 on other levels

  if (ilevel == 0) foriginal_saved[0] = foriginal_saved[1] = foriginal_saved[2] = 0.0;

  if (ilevel == ilevel_respa) {
    post_force(vflag);
    foriginal[0] += foriginal_saved[0];
    foriginal[1] += foriginal_saved[1];
    foriginal[2] += foriginal_saved[2];
  } else {
    if (region) region->prematch();

    double **x = atom->x;
    double **f = atom->f;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        foriginal_saved[0] += f[i][0];
        foriginal_saved[1] += f[i][1];
        foriginal_saved[2] += f[i][2];
        if (xstyle) f[i][0] = 0.0;
        if (ystyle) f[i][1] = 0.0;
        if (zstyle) f[i][2] = 0.0;
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixSetForce::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total force on fix group before force was changed
------------------------------------------------------------------------- */

double FixSetForce::compute_vector(int n)
{
  // only sum across procs one time

  if (force_flag == 0) {
    MPI_Allreduce(foriginal, foriginal_all, 3, MPI_DOUBLE, MPI_SUM, world);
    force_flag = 1;
  }
  return foriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixSetForce::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom * 3 * sizeof(double);
  return bytes;
}
