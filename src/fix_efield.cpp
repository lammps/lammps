/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Christina Payne (Vanderbilt U)
------------------------------------------------------------------------- */

#include "string.h"
#include "stdlib.h"
#include "fix_efield.h"
#include "atom.h"
#include "update.h"
#include "force.h"
#include "respa.h"
#include "error.h"
#include "math.h"

/* ---------------------------------------------------------------------- */

FixEfield::FixEfield(int narg, char **arg) : Fix(narg, arg)
{
  if (narg != 6) error->all("Illegal fix efield command");

  double factor = force->qe2f;
  ex = factor * atof(arg[3]);
  ey = factor * atof(arg[4]);
  ez = factor * atof(arg[5]);
}

/* ---------------------------------------------------------------------- */

int FixEfield::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEfield::init()
{
  // require an atom style with charge defined

  if (atom->charge_allow == 0)
    error->all("Must use charged atom style with fix efield");

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixEfield::setup()
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(1,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ----------------------------------------------------------------------
   apply F = qE
------------------------------------------------------------------------- */

void FixEfield::post_force(int vflag)
{
  double **f = atom->f;
  double *q = atom->q;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0] += q[i]*ex;
      f[i][1] += q[i]*ey;
      f[i][2] += q[i]*ez;
    }
}

/* ---------------------------------------------------------------------- */

void FixEfield::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}
