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

#include "fix_addtorque_atom.h"

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
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

FixAddTorqueAtom::FixAddTorqueAtom(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), xstr(nullptr), ystr(nullptr), zstr(nullptr),
    idregion(nullptr), region(nullptr), storque(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix addtorque/atom", error);

  dynamic_group_allow = 1;
  vector_flag = 1;
  size_vector = 3;
  global_freq = 1;
  extvector = 1;
  respa_level_support = 1;
  ilevel_respa = 0;

  xstyle = ystyle = zstyle = NONE;
  xvar = yvar = zvar = -1;
  varflag = NONE;

  if (utils::strmatch(arg[3], "^v_")) {
    xstr = utils::strdup(arg[3] + 2);
  } else {
    xvalue = utils::numeric(FLERR, arg[3], false, lmp);
    xstyle = CONSTANT;
  }
  if (utils::strmatch(arg[4], "^v_")) {
    ystr = utils::strdup(arg[4] + 2);
  } else {
    yvalue = utils::numeric(FLERR, arg[4], false, lmp);
    ystyle = CONSTANT;
  }
  if (utils::strmatch(arg[5], "^v_")) {
    zstr = utils::strdup(arg[5] + 2);
  } else {
    zvalue = utils::numeric(FLERR, arg[5], false, lmp);
    zstyle = CONSTANT;
  }

  // optional args

  nevery = 1;
  int iarg = 6;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "every") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix addtorque/atom every", error);
      nevery = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (nevery <= 0) error->all(FLERR, "Invalid fix addtorque/atom every argument: {}", nevery);
      iarg += 2;
    } else if (strcmp(arg[iarg], "region") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix addtorque/atom region", error);
      region = domain->get_region_by_id(arg[iarg + 1]);
      if (!region) error->all(FLERR, "Region {} for fix addtorque/atom does not exist", arg[iarg + 1]);
      delete[] idregion;
      idregion = utils::strdup(arg[iarg + 1]);
      iarg += 2;
    } else
      error->all(FLERR, "Unknown fix addtorque/atom keyword: {}", arg[iarg]);
  }

  torque_flag = 0;
  toriginal[0] = toriginal[1] = toriginal[2] = 0.0;

  maxatom = 1;
  memory->create(storque, maxatom, 3, "addtorque/atom:storque");
}

/* ---------------------------------------------------------------------- */

FixAddTorqueAtom::~FixAddTorqueAtom()
{
  delete[] xstr;
  delete[] ystr;
  delete[] zstr;
  delete[] idregion;
  memory->destroy(storque);
}

/* ---------------------------------------------------------------------- */

int FixAddTorqueAtom::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAddTorqueAtom::init()
{
  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR, "Variable {} for fix addtorque/atom does not exist", xstr);
    if (input->variable->equalstyle(xvar))
      xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar))
      xstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix addtorque/atom is invalid style", xstr);
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR, "Variable {} for fix addtorque/atom does not exist", ystr);
    if (input->variable->equalstyle(yvar))
      ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar))
      ystyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix addtorque/atom is invalid style", ystr);
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR, "Variable {} for fix addtorque/atom does not exist", zstr);
    if (input->variable->equalstyle(zvar))
      zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar))
      zstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix addtorque/atom is invalid style", zstr);
  }

  // set index and check validity of region

  if (idregion) {
    region = domain->get_region_by_id(idregion);
    if (!region) error->all(FLERR, "Region {} for fix addtorque/atom does not exist", idregion);
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else
    varflag = CONSTANT;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, ilevel_respa);
  }

  if ((modify->check_rigid_group_overlap(groupbit)) && (comm->me == 0))
    error->warning(FLERR,"Adding torques to atoms in rigid bodies with fix addtorque/atom may not work as expected");
}

/* ---------------------------------------------------------------------- */

void FixAddTorqueAtom::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag, ilevel_respa, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixAddTorqueAtom::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddTorqueAtom::post_force(int /*vflag*/)
{
  double **x = atom->x;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (update->ntimestep % nevery) return;

  // update region if necessary

  if (region) region->prematch();

  // reallocate storque array if necessary

  if ((varflag == ATOM) && atom->nmax > maxatom) {
    maxatom = atom->nmax;
    memory->destroy(storque);
    memory->create(storque, maxatom, 3, "addtorque/atom:storque");
  }

  // toriginal[012] = torque on atoms before extra torque added

  toriginal[0] = toriginal[1] = toriginal[2] = 0.0;
  torque_flag = 0;

  // constant torque

  if (varflag == CONSTANT) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        toriginal[0] += torque[i][0];
        toriginal[1] += torque[i][1];
        toriginal[2] += torque[i][2];
        torque[i][0] += xvalue;
        torque[i][1] += yvalue;
        torque[i][2] += zvalue;
      }

    // variable torque, wrap with clear/add
    // wrap with clear/add

  } else {
    modify->clearstep_compute();

    if (xstyle == EQUAL)
      xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar, igroup, &storque[0][0], 3, 0);
    if (ystyle == EQUAL)
      yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar, igroup, &storque[0][1], 3, 0);
    if (zstyle == EQUAL)
      zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar, igroup, &storque[0][2], 3, 0);

    modify->addstep_compute(update->ntimestep + 1);

    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
        if (xstyle == ATOM) xvalue = storque[i][0];
        if (ystyle == ATOM) yvalue = storque[i][1];
        if (zstyle == ATOM) zvalue = storque[i][2];

        toriginal[0] += torque[i][0];
        toriginal[1] += torque[i][1];
        toriginal[2] += torque[i][2];

        if (xstyle) torque[i][0] += xvalue;
        if (ystyle) torque[i][1] += yvalue;
        if (zstyle) torque[i][2] += zvalue;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixAddTorqueAtom::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAddTorqueAtom::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   return components of total torque on fix group before torque was changed
------------------------------------------------------------------------- */

double FixAddTorqueAtom::compute_vector(int n)
{
  // only sum across procs one time

  if (torque_flag == 0) {
    MPI_Allreduce(toriginal, toriginal_all, 3, MPI_DOUBLE, MPI_SUM, world);
    torque_flag = 1;
  }
  return toriginal_all[n];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixAddTorqueAtom::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = maxatom * 3 * sizeof(double);
  return bytes;
}
