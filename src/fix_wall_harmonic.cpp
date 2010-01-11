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
#include "fix_wall_harmonic.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallHarmonic::FixWallHarmonic(LAMMPS *lmp, int narg, char **arg) : 
  FixWall(lmp, narg, arg) {}

/* ----------------------------------------------------------------------
   interaction of all particles in group with all 6 walls (if defined)
   error if any particle is on or behind wall
------------------------------------------------------------------------- */

void FixWallHarmonic::wall_particle(int m, double coord)
{
  double delta,dr,fwall;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int dim = m/2;
  int side = m % 2;
  if (side == 0) side = -1;

  int onflag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side < 0) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta >= cutoff[m]) continue;
      if (delta <= 0.0) {
	onflag = 1;
	continue;
      }
      dr = cutoff[m]-delta;
      fwall = side * 2.0*epsilon[m]*dr;
      f[i][dim] -= fwall;
      ewall[0] += epsilon[m]*dr*dr;
      ewall[m+1] += fwall;
    }

  if (onflag) error->one("Particle on or inside fix wall surface");
}
