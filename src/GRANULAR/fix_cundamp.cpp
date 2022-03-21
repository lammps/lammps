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

#include "fix_cundamp.h"

#include "atom.h"
#include "error.h"
#include "respa.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixCundamp::FixCundamp(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  gamma_lin(nullptr),gamma_ang(nullptr)
{
  dynamic_group_allow = 1;

  if (!atom->sphere_flag)
    error->all(FLERR,"Fix cundamp requires atom style sphere");

  if (narg < 5) error->all(FLERR,"Illegal fix cundamp command");

  double gamma_lin_one = utils::numeric(FLERR,arg[3],false,lmp);
  double gamma_ang_one = utils::numeric(FLERR,arg[4],false,lmp);
  gamma_lin = new double[atom->ntypes+1];
  gamma_ang = new double[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++) {
    gamma_lin[i] = gamma_lin_one;
    gamma_ang[i] = gamma_ang_one;
  }

  // optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix cundamp command");
      int itype = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      double scale = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (itype <= 0 || itype > atom->ntypes)
        error->all(FLERR,"Illegal fix cundamp command");
      gamma_lin[itype] = gamma_lin_one * scale;
      gamma_ang[itype] = gamma_ang_one * scale;
      iarg += 3;
    } else error->all(FLERR,"Illegal fix cundamp command");
  }

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixCundamp::~FixCundamp()
{
  delete [] gamma_lin;
  delete [] gamma_ang;
}

/* ---------------------------------------------------------------------- */

int FixCundamp::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCundamp::init()
{
  int max_respa = 0;

  if (utils::strmatch(update->integrate_style,"^respa")) {
    ilevel_respa = max_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,max_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixCundamp::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixCundamp::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCundamp::post_force(int /*vflag*/)
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

  double gamma_l,gamma_a;
  int signf0,signf1,signf2;
  int signt0,signt1,signt2;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      gamma_l = gamma_lin[type[i]];
      gamma_a = gamma_ang[type[i]];

      signf0 = (f[i][0]*v[i][0] >= 0.0) ? 1.0 : -1.0;
      signf1 = (f[i][1]*v[i][1] >= 0.0) ? 1.0 : -1.0;
      signf2 = (f[i][2]*v[i][2] >= 0.0) ? 1.0 : -1.0;
      f[i][0] *= 1.0 - gamma_l*signf0;
      f[i][1] *= 1.0 - gamma_l*signf1;
      f[i][2] *= 1.0 - gamma_l*signf2;

      signt0 = (torque[i][0]*omega[i][0] >= 0.0) ? 1.0 : -1.0;
      signt1 = (torque[i][1]*omega[i][1] >= 0.0) ? 1.0 : -1.0;
      signt2 = (torque[i][2]*omega[i][2] >= 0.0) ? 1.0 : -1.0;
      torque[i][0] *= 1.0 - gamma_a*signt0;
      torque[i][1] *= 1.0 - gamma_a*signt1;
      torque[i][2] *= 1.0 - gamma_a*signt2;
    }
}

/* ---------------------------------------------------------------------- */

void FixCundamp::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCundamp::min_post_force(int vflag)
{
  post_force(vflag);
}
