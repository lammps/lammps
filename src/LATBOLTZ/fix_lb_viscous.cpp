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

/* ----------------------------------------------------------------------
   Contributing authors: Frances Mackay, Santtu Ollila, Colin Denniston (UWO)
------------------------------------------------------------------------- */

#include "fix_lb_viscous.h"

#include "atom.h"
#include "error.h"
#include "fix_lb_fluid.h"
#include "group.h"
#include "modify.h"
#include "respa.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLbViscous::FixLbViscous(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), fix_lb_fluid(nullptr)
{
  if (narg < 3) error->all(FLERR, "Illegal fix lb/viscous command");

  int groupbit_lb_fluid = 0;

  auto ifix = modify->get_fix_by_style("lb/fluid");
  if (ifix.size() > 0) {
    fix_lb_fluid = (FixLbFluid *) ifix[0];
    groupbit_lb_fluid = group->bitmask[ifix[0]->igroup];
  } else {
    error->all(FLERR, "the lb/fluid fix must also be used if using the lb/viscous fix");
  }

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  for (int j = 0; j < nlocal; j++) {
    if ((mask[j] & groupbit) && !(mask[j] & groupbit_lb_fluid))
      error->one(FLERR, "Atoms must be in the fix lb/fluid group to apply a fluid force to them");
  }
}

/* ---------------------------------------------------------------------- */

int FixLbViscous::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLbViscous::init()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixLbViscous::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa - 1);
    post_force_respa(vflag, nlevels_respa - 1, 0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa - 1);
  }
}

/* ---------------------------------------------------------------------- */

void FixLbViscous::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLbViscous::post_force(int /*vflag*/)
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type

  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;
  double *massp = fix_lb_fluid->massp;
  double massfactor;
  double **hydroF = fix_lb_fluid->hydroF;

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massfactor = rmass[i] / (rmass[i] + massp[i]);

        f[i][0] = hydroF[i][0] + massfactor * f[i][0];
        f[i][1] = hydroF[i][1] + massfactor * f[i][1];
        f[i][2] = hydroF[i][2] + massfactor * f[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        massfactor = mass[type[i]] / (mass[type[i]] + massp[i]);

        f[i][0] = hydroF[i][0] + massfactor * f[i][0];
        f[i][1] = hydroF[i][1] + massfactor * f[i][1];
        f[i][2] = hydroF[i][2] + massfactor * f[i][2];
      }
  }
}

/* ---------------------------------------------------------------------- */

void FixLbViscous::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLbViscous::min_post_force(int vflag)
{
  post_force(vflag);
}
