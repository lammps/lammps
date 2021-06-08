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

/* ----------------------------------------------------------------------
   Contributing author: Zheng GONG (ENS de Lyon, z.gong@outlook.com)
------------------------------------------------------------------------- */

#include "fix_accelerate_cos.h"

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixAccelerateCos::FixAccelerateCos(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg) {
  if (narg < 4) error->all(FLERR, "Illegal fix accelerate/cos command");
  acceleration = utils::numeric(FLERR, arg[3],false,lmp);
  if (domain->dimension == 2)
    error->all(FLERR,"Fix accelerate/cos cannot be used with 2d systems");
}

/* ---------------------------------------------------------------------- */

int FixAccelerateCos::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAccelerateCos::setup(int vflag) {
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAccelerateCos::post_force(int /* vflag */) {
  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double massone, force_x, acc_x;
  double zlo = domain->boxlo[2];
  double zhi = domain->boxhi[2];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      acc_x = acceleration *
              cos(MathConst::MY_2PI * (x[i][2] - zlo) / (zhi - zlo));
      force_x = acc_x * massone * force->mvv2e;

      f[i][0] += force_x;
    }
}
