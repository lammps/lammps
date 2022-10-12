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

#include "fix_wall_morse.h"
#include "atom.h"
#include "error.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallMorse::FixWallMorse(LAMMPS *lmp, int narg, char **arg) : FixWall(lmp, narg, arg)
{
  dynamic_group_allow = 1;
}

/* ---------------------------------------------------------------------- */

void FixWallMorse::precompute(int m)
{
  coeff1[m] = 2.0 * epsilon[m] * alpha[m];
  const double alpha_dr = -alpha[m] * (cutoff[m] - sigma[m]);
  offset[m] = epsilon[m] * (exp(2.0 * alpha_dr) - 2.0 * exp(alpha_dr));
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with a wall
   m = index of wall coeffs
   which = xlo,xhi,ylo,yhi,zlo,zhi
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallMorse::wall_particle(int m, int which, double coord)
{
  double delta, fwall;
  double vn;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      if (side < 0)
        delta = x[i][dim] - coord;
      else
        delta = coord - x[i][dim];
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
        onflag = 1;
        continue;
      }
      double dr = delta - sigma[m];
      double dexp = exp(-alpha[m] * dr);
      fwall = side * coeff1[m] * (dexp * dexp - dexp) / delta;
      ewall[0] += epsilon[m] * (dexp * dexp - 2.0 * dexp) - offset[m];
      f[i][dim] -= fwall;
      ewall[m + 1] += fwall;

      if (evflag) {
        if (side < 0)
          vn = -fwall * delta;
        else
          vn = fwall * delta;
        v_tally(dim, i, vn);
      }
    }
  }

  if (onflag) error->one(FLERR, "Particle on or inside fix wall surface");
}
