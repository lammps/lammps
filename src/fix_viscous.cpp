/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_viscous.h"
#include <cstring>
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixViscous::FixViscous(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg),
  gamma(NULL)
{
  dynamic_group_allow = 1;

  if (narg < 4) error->all(FLERR,"Illegal fix viscous command");

  double gamma_one = force->numeric(FLERR,arg[3]);
  gamma = new double[atom->ntypes+1];
  for (int i = 1; i <= atom->ntypes; i++) gamma[i] = gamma_one;

  // optional args

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix viscous command");
      int itype = force->inumeric(FLERR,arg[iarg+1]);
      double scale = force->numeric(FLERR,arg[iarg+2]);
      if (itype <= 0 || itype > atom->ntypes)
        error->all(FLERR,"Illegal fix viscous command");
      gamma[itype] = gamma_one * scale;
      iarg += 3;
    } else error->all(FLERR,"Illegal fix viscous command");
  }

  respa_level_support = 1;
  ilevel_respa = 0;
}

/* ---------------------------------------------------------------------- */

FixViscous::~FixViscous()
{
  delete [] gamma;
}

/* ---------------------------------------------------------------------- */

int FixViscous::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixViscous::init()
{
  int max_respa = 0;

  if (strstr(update->integrate_style,"respa")) {
    ilevel_respa = max_respa = ((Respa *) update->integrate)->nlevels-1;
    if (respa_level >= 0) ilevel_respa = MIN(respa_level,max_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixViscous::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(ilevel_respa);
    post_force_respa(vflag,ilevel_respa,0);
    ((Respa *) update->integrate)->copy_f_flevel(ilevel_respa);
  }
}

/* ---------------------------------------------------------------------- */

void FixViscous::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscous::post_force(int /*vflag*/)
{
  // apply drag force to atoms in group
  // direction is opposed to velocity vector
  // magnitude depends on atom type

  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double drag;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      drag = gamma[type[i]];
      f[i][0] -= drag*v[i][0];
      f[i][1] -= drag*v[i][1];
      f[i][2] -= drag*v[i][2];
    }
}

/* ---------------------------------------------------------------------- */

void FixViscous::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == ilevel_respa) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixViscous::min_post_force(int vflag)
{
  post_force(vflag);
}
