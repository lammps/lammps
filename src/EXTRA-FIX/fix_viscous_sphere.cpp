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

#include "fix_viscous_sphere.h"

#include "atom.h"
#include "error.h"
#include "input.h"
#include "memory.h"
#include "respa.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

// type of scaling

enum { NONE, TYPE, VARIABLE };

/* ---------------------------------------------------------------------- */

FixViscousSphere::FixViscousSphere(LAMMPS *_lmp, int narg, char **arg) :
    Fix(_lmp, narg, arg), scalegamma(nullptr), scaleval(nullptr), scalevarid(nullptr),
    scalestyle(NONE), scalevar(-1)
{
  dynamic_group_allow = 1;

  if (!atom->sphere_flag) error->all(FLERR, "Fix viscous/sphere requires atom style sphere");

  if (narg < 4) error->all(FLERR, "Illegal fix viscous/sphere command");

  gamma = utils::numeric(FLERR, arg[3], false, lmp);

  // optional args

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "scale") == 0) {
      if (iarg + 2 > narg) error->all(FLERR, "Illegal fix viscous/sphere command");
      if (utils::strmatch(arg[iarg + 1], "^v_")) {
        if (scalestyle != NONE) error->all(FLERR, "Must use only one style of scaling");
        scalevarid = utils::strdup(arg[iarg + 1] + 2);
        int ivariable = input->variable->find(scalevarid);
        if (ivariable < 0)
          error->all(FLERR, "Variable name {} for fix viscous/sphere does not exist", scalevarid);
        if (input->variable->atomstyle(ivariable) == 0)
          error->all(FLERR, "Fix viscous/scale variable {} is not atom-style variable", scalevarid);
        scalestyle = VARIABLE;
        memory->destroy(scaleval);
        memory->create(scaleval, atom->nmax, "fix_viscous/sphere:scaleval");
        iarg += 2;
      } else {
        if (scalestyle == VARIABLE) error->all(FLERR, "Must use only one style of scaling");
        if (iarg + 3 > narg) error->all(FLERR, "Illegal fix viscous/sphere command");
        if (!scalegamma) {
          scalegamma = new double[atom->ntypes + 1];
          for (int i = 1; i <= atom->ntypes; i++) scalegamma[i] = 1.0;
        }
        int itype = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
        double scale = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
        if ((itype < 1) || (itype > atom->ntypes))
          error->all(FLERR, "Atom type {} out of range for fix viscous/sphere command:", itype);
        scalegamma[itype] = scale;
        scalestyle = TYPE;
        iarg += 3;
      }
    } else
      error->all(FLERR, "Illegal fix viscous/sphere command");
  }

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixViscousSphere::~FixViscousSphere()
{
  memory->destroy(scaleval);
  delete[] scalegamma;
  delete[] scalevarid;
}

/* ---------------------------------------------------------------------- */

int FixViscousSphere::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixViscousSphere::init()
{
  int max_respa = 0;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = max_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, max_respa);
  }

  if (scalestyle == VARIABLE) {
    int ivariable = input->variable->find(scalevarid);
    if (ivariable < 0)
      error->all(FLERR, "Variable name {} for fix viscous/sphere does not exist", scalevarid);
    if (input->variable->atomstyle(ivariable) == 0)
      error->all(FLERR, "Fix viscous/sphere variable {} is not atom-style variable", scalevarid);
    scalevar = ivariable;
  }
}

/* ---------------------------------------------------------------------- */

void FixViscousSphere::setup(int vflag)
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

void FixViscousSphere::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousSphere::post_force(int /*vflag*/)
{
  // apply drag torque to finite-size atoms in group
  // direction is opposed to angular velocity vector
  // magnitude depends on atom type

  double **omega = atom->omega;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double drag;
  if (scalestyle == VARIABLE) {
    memory->grow(scaleval, atom->nmax, "fix_viscous/sphere:scaleval");
    input->variable->compute_atom(scalevar, igroup, scaleval, 1, 0);
  }

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (scalestyle == TYPE) {
        drag = gamma * scalegamma[type[i]];
      } else if (scalestyle == VARIABLE) {
        drag = gamma * scaleval[i];
      } else {    // scalestyle == NONE
        drag = gamma;
      }
      torque[i][0] -= drag * omega[i][0];
      torque[i][1] -= drag * omega[i][1];
      torque[i][2] -= drag * omega[i][2];
    }
}

/* ---------------------------------------------------------------------- */

void FixViscousSphere::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousSphere::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixViscousSphere::memory_usage()
{
  if (scalestyle == VARIABLE)
    return (double) sizeof(double) * atom->nmax;
  else if (scalestyle == TYPE)
    return (double) sizeof(double) * atom->ntypes;
  else
    return 0.0;
}
