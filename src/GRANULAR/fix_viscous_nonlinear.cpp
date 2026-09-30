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

#include "fix_viscous_nonlinear.h"

#include "atom.h"
#include "error.h"
#include "math_const.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

/* ----------------------------------------------------------------------
   Nonlinear (Schiller-Naumann) fluid drag on finite-size spherical particles:

     F = -1/2 C_d rho_f (pi r^2) |v_rel| v_rel,   v_rel = v_particle - v_fluid
     C_d = 24/Re (1 + 0.15 Re^0.687),   Re = rho_f |v_rel| (2 r) / mu_f

   In the low-Reynolds limit this reduces to Stokes drag 6 pi mu_f r v_rel.
------------------------------------------------------------------------- */

FixViscousNonlinear::FixViscousNonlinear(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg)
{
  dynamic_group_allow = 1;

  if (narg < 5) error->all(FLERR, "Illegal fix viscous/nonlinear command");

  rho_fluid = utils::numeric(FLERR, arg[3], false, lmp);
  mu_fluid = utils::numeric(FLERR, arg[4], false, lmp);
  if (rho_fluid <= 0.0 || mu_fluid <= 0.0)
    error->all(FLERR, "Fix viscous/nonlinear fluid density and viscosity must be > 0");

  v_fluid[0] = v_fluid[1] = v_fluid[2] = 0.0;

  // optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "velocity") == 0) {
      if (iarg + 4 > narg) utils::missing_cmd_args(FLERR, "fix viscous/nonlinear velocity", error);
      v_fluid[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      v_fluid[1] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      v_fluid[2] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else
      error->all(FLERR, "Illegal fix viscous/nonlinear command");
  }

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixViscousNonlinear::~FixViscousNonlinear()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

int FixViscousNonlinear::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixViscousNonlinear::init()
{
  if (!atom->radius_flag)
    error->all(FLERR, "Fix viscous/nonlinear requires atom attribute radius");

  int max_respa = 0;

  if (utils::strmatch(update->integrate_style, "^respa")) {
    ilevel_respa = max_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels - 1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level, max_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixViscousNonlinear::setup(int vflag)
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

void FixViscousNonlinear::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousNonlinear::post_force(int /*vflag*/)
{
  // apply Schiller-Naumann drag relative to the (uniform) fluid velocity

  double **v = atom->v;
  double **f = atom->f;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double vrel[3];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      vrel[0] = v[i][0] - v_fluid[0];
      vrel[1] = v[i][1] - v_fluid[1];
      vrel[2] = v[i][2] - v_fluid[2];
      const double vmag = sqrt(vrel[0] * vrel[0] + vrel[1] * vrel[1] + vrel[2] * vrel[2]);
      if (vmag == 0.0) continue;

      const double r = radius[i];
      const double re = rho_fluid * vmag * (2.0 * r) / mu_fluid;
      const double cd = (24.0 / re) * (1.0 + 0.15 * pow(re, 0.687));

      // F = -1/2 Cd rho_g (pi r^2) |v_rel| v_rel
      const double pref = 0.5 * cd * rho_fluid * MY_PI * r * r * vmag;
      f[i][0] -= pref * vrel[0];
      f[i][1] -= pref * vrel[1];
      f[i][2] -= pref * vrel[2];
    }
}

/* ---------------------------------------------------------------------- */

void FixViscousNonlinear::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscousNonlinear::min_post_force(int vflag)
{
  post_force(vflag);
}
