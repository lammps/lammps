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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_lineforce.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLineForce::FixLineForce(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 6) error->all(FLERR,"Illegal fix lineforce command");
  xdir = force->numeric(FLERR,arg[3]);
  ydir = force->numeric(FLERR,arg[4]);
  zdir = force->numeric(FLERR,arg[5]);

  double len = sqrt(xdir*xdir + ydir*ydir + zdir*zdir);
  if (len == 0.0) error->all(FLERR,"Illegal fix lineforce command");

  xdir /= len;
  ydir /= len;
  zdir /= len;
}

/* ---------------------------------------------------------------------- */

int FixLineForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLineForce::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLineForce::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLineForce::post_force(int vflag)
{
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dot;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dot = f[i][0]*xdir + f[i][1]*ydir + f[i][2]*zdir;
      f[i][0] = dot * xdir;
      f[i][1] = dot * ydir;
      f[i][2] = dot * zdir;
    }
}

/* ---------------------------------------------------------------------- */

void FixLineForce::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLineForce::min_post_force(int vflag)
{
  post_force(vflag);
}
