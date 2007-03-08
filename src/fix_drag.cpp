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
#include "stdlib.h"
#include "string.h"
#include "fix_drag.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "domain.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixDrag::FixDrag(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all("Illegal fix drag command");

  xflag = yflag = zflag = 1;

  if (strcmp(arg[3],"NULL") == 0) xflag = 0;
  else xc = atof(arg[3]);
  if (strcmp(arg[4],"NULL") == 0) yflag = 0;
  else yc = atof(arg[4]);
  if (strcmp(arg[5],"NULL") == 0) zflag = 0;
  else zc = atof(arg[5]);

  f_mag = atof(arg[6]);
  delta = atof(arg[7]);
}

/* ---------------------------------------------------------------------- */

int FixDrag::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixDrag::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixDrag::setup()
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(1,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixDrag::post_force(int vflag)
{
  // apply drag force to atoms in group of magnitude f_mag
  // apply in direction (r-r0) if atom is further than delta away

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  
  double dx,dy,dz,r,prefactor;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dx = x[i][0] - xc;
      dy = x[i][1] - yc;
      dz = x[i][2] - zc;
      if (!xflag) dx = 0.0;
      if (!yflag) dy = 0.0;
      if (!zflag) dz = 0.0;
      domain->minimum_image(dx,dy,dz);
      r = sqrt(dx*dx + dy*dy + dz*dz);
      if (r > delta) {
	prefactor = f_mag/r;
	f[i][0] -= prefactor*dx;
	f[i][1] -= prefactor*dy;
	f[i][2] -= prefactor*dz;
      }
    }
}

/* ---------------------------------------------------------------------- */

void FixDrag::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}
