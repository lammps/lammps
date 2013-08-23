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

/* ----------------------------------------------------------------------
   Contributing author: Jonathan Lee (Sandia)
------------------------------------------------------------------------- */

#include "math.h"
#include "fix_wall_lj1043.h"
#include "atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallLJ1043::FixWallLJ1043(LAMMPS *lmp, int narg, char **arg) : 
  FixWall(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixWallLJ1043::precompute(int m)
{
  coeff1[m] = 2.0/5.0 * epsilon[m] * pow(sigma[m],10.0);
  coeff2[m] = epsilon[m] * pow(sigma[m],4.0);
  coeff3[m] = pow(2.0,1/2.0) / 3 * epsilon[m] * pow(sigma[m],3.0);
  coeff4[m] = 0.61 / pow(2.0,1/2.0) * sigma[m];
  coeff5[m] = coeff1[m] * 10.0;
  coeff6[m] = coeff2[m] * 4.0;
  coeff7[m] = coeff3[m] * 3.0;

  double rinv = 1.0/cutoff[m];
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  offset[m] = coeff1[m]*r4inv*r4inv*r2inv - coeff2[m]*r4inv - 
	coeff3[m]*pow(cutoff[m]+coeff4[m],-3.0);
}

/* ---------------------------------------------------------------------- */

void FixWallLJ1043::wall_particle(int m, int which, double coord)
{
  double delta,rinv,r2inv,r4inv,r10inv,fwall;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = which / 2;
  int side = which % 2;
  if (side == 0) side = -1;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta <= 0.0) continue;
      if (delta > cutoff[m]) continue;
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r10inv = r4inv*r4inv*r2inv;

      fwall = side * (coeff5[m]*r10inv*rinv - coeff6[m]*r4inv*rinv - 
	coeff7[m]*pow(delta+coeff4[m],-4.0));
      f[i][dim] -= fwall;
      ewall[0] += coeff1[m]*r10inv - coeff2[m]*r4inv - 
	coeff3[m]*pow(delta+coeff4[m],-3.0) - offset[m];
      ewall[m+1] += fwall;
    }
}
