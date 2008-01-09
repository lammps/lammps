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

#include "string.h"
#include "fix_freeze.h"
#include "atom.h"
#include "modify.h"
#include "comm.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixFreeze::FixFreeze(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal fix freeze command");

  if (!atom->torque_flag)
    error->all("Fix freeze requires atom attribute torque");
}

/* ---------------------------------------------------------------------- */

int FixFreeze::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFreeze::init()
{
  // error if more than one freeze fix
  // because accessed by pair style granular and fix gran/diag

  int count = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) count++;
  if (count > 1) error->all("More than one fix freeze");
}

/* ---------------------------------------------------------------------- */

void FixFreeze::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFreeze::post_force(int vflag)
{
  double **f = atom->f;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      f[i][0] = 0.0;
      f[i][1] = 0.0;
      f[i][2] = 0.0;
      torque[i][0] = 0.0;
      torque[i][1] = 0.0;
      torque[i][2] = 0.0;
    }
}
