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

#include "fix_damping_cundall.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// type of scaling

enum { NONE, TYPE, VARIABLE };

/* ---------------------------------------------------------------------- */

FixDampingCundall::FixDampingCundall(LAMMPS *_lmp, int narg, char **arg) :
    Fix(_lmp, narg, arg), scalegamma(nullptr), scaleval(nullptr), scalevarid(nullptr),
    scalestyle(NONE), scalevar(-1)
{
  dynamic_group_allow = 1;

  if (!atom->omega_flag) error->all(FLERR, "Fix damping/cundall requires atom attribute omega");

  if (narg < 5) utils::missing_cmd_args(FLERR, "fix damping/cundall", error);

  gamma_lin = utils::numeric(FLERR, arg[3], false, lmp);
  gamma_ang = utils::numeric(FLERR, arg[4], false, lmp);

  // optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "scale") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix damping/cundall command");
      if (utils::strmatch(arg[iarg + 1], "^v_")) {
        if (scalestyle != NONE) error->all(FLERR, "Must use only one style of scaling");
        scalevarid = utils::strdup(arg[iarg + 1] + 2);
        int ivariable = input->variable->find(scalevarid);
        if (ivariable < 0)
          error->all(FLERR, "Variable name {} for fix damping/cundall does not exist", scalevarid);
        if (input->variable->atomstyle(ivariable) == 0)
          error->all(FLERR, "Fix viscous/scale variable {} is not atom-style variable", scalevarid);
        scalestyle = VARIABLE;
        memory->destroy(scaleval);
        memory->create(scaleval, atom->nmax, "fix_damping/cundall:scaleval");
        iarg += 2;
      } else {
        if (scalestyle == VARIABLE) error->all(FLERR, "Must use only one style of scaling");
        if (iarg + 3 > narg) error->all(FLERR, "Illegal fix damping/cundall command");
        if (!scalegamma) {
          scalegamma = new double[atom->ntypes + 1];
          for (int i = 1; i <= atom->ntypes; i++) scalegamma[i] = 1.0;
        }
        int itype = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
        double scale = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        if ((itype < 1) || (itype > atom->ntypes))
          error->all(FLERR, "Atom type {} out of range for fix damping/cundall command:", itype);
        scalegamma[itype] = scale;
        scalestyle = TYPE;
        iarg += 3;
      }
    } else
      error->all(FLERR, "Illegal fix damping/cundall command");
  }

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixDampingCundall::~FixDampingCundall()
{
  memory->destroy(scaleval);
  delete[] scalegamma;
  delete[] scalevarid;
}

/* ---------------------------------------------------------------------- */

int FixDampingCundall::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDampingCundall::init()
{
  int max_respa = 0;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = max_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, max_respa);
  }

  bool fflag = false;
  for (auto &ifix : modify->get_fix_list()) {
    if (fflag && (comm->me == 0) && (ifix->setmask() & POST_FORCE))
      error->warning(FLERR, "Fix {} alters forces after fix damping/cundall", ifix->id);
    if (ifix == this) fflag = true;
  }

  if (scalestyle == VARIABLE) {
    int ivariable = input->variable->find(scalevarid);
    if (ivariable < 0)
      error->all(FLERR, "Variable name {} for fix damping/cundall does not exist", scalevarid);
    if (input->variable->atomstyle(ivariable) == 0)
      error->all(FLERR, "Fix damping/cundall variable {} is not atom-style variable", scalevarid);
    scalevar = ivariable;
  }
}

/* ---------------------------------------------------------------------- */

void FixDampingCundall::setup(int vflag)
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

void FixDampingCundall::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixDampingCundall::post_force(int /*vflag*/)
{
  // apply damping force/torque to finite-size atoms in group
  // add a fraction of the current force/torque if work is negative
  // subtract a fraction of the current force/torque if work is positive
  // applied over each component independently (non-objective)
  // magnitude depends on atom type

  // see, e.g. Yade-DEM implementation of NewtonIntegrator::cundallDamp1st()
  // gitlab.com/yade-dev/trunk/-/blob/master/pkg/dem/NewtonIntegrator.cpp

  double **v = atom->v;
  double **omega = atom->omega;
  double **f = atom->f;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double gamma_l, gamma_a;
  int signf0, signf1, signf2;
  int signt0, signt1, signt2;

  if (scalestyle == VARIABLE) {
    memory->grow(scaleval, atom->nmax, "fix_damping/cundall:scaleval");
    input->variable->compute_atom(scalevar, igroup, scaleval, 1, 0);
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (scalestyle == TYPE) {
        gamma_l = gamma_lin * scalegamma[type[i]];
        gamma_a = gamma_ang * scalegamma[type[i]];
      } else if (scalestyle == VARIABLE) {
        gamma_l = gamma_lin * scaleval[i];
        gamma_a = gamma_ang * scaleval[i];
      } else {    // scalestyle NONE
        gamma_l = gamma_lin;
        gamma_a = gamma_ang;
      }

      signf0 = (f[i][0] * v[i][0] >= 0.0) ? 1.0 : -1.0;
      signf1 = (f[i][1] * v[i][1] >= 0.0) ? 1.0 : -1.0;
      signf2 = (f[i][2] * v[i][2] >= 0.0) ? 1.0 : -1.0;
      f[i][0] *= 1.0 - gamma_l * signf0;
      f[i][1] *= 1.0 - gamma_l * signf1;
      f[i][2] *= 1.0 - gamma_l * signf2;

      signt0 = (torque[i][0] * omega[i][0] >= 0.0) ? 1.0 : -1.0;
      signt1 = (torque[i][1] * omega[i][1] >= 0.0) ? 1.0 : -1.0;
      signt2 = (torque[i][2] * omega[i][2] >= 0.0) ? 1.0 : -1.0;
      torque[i][0] *= 1.0 - gamma_a * signt0;
      torque[i][1] *= 1.0 - gamma_a * signt1;
      torque[i][2] *= 1.0 - gamma_a * signt2;
    }
}

/* ---------------------------------------------------------------------- */

void FixDampingCundall::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixDampingCundall::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixDampingCundall::memory_usage()
{
  if (scalestyle == VARIABLE)
    return (double) sizeof(double) * atom->nmax;
  else if (scalestyle == TYPE)
    return 2.0 * sizeof(double) * atom->ntypes;
  else
    return 0.0;
}
