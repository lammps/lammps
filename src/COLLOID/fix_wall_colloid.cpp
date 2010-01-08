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
   Contributing authors: Jeremy Lechman (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "fix_wall_colloid.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallColloid::FixWallColloid(LAMMPS *lmp, int narg, char **arg) : 
  FixWall(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixWallColloid::init()
{
  if (!atom->avec->shape_type)
    error->all("Fix wall/colloid requires atom attribute shape");
  if (atom->radius_flag)
    error->all("Fix wall/colloid cannot be used with atom attribute diameter");

  // insure all particle shapes are spherical
  // can be polydisperse

  for (int i = 1; i <= atom->ntypes; i++)
    if ((atom->shape[i][0] != atom->shape[i][1]) || 
	(atom->shape[i][0] != atom->shape[i][2]) ||
	(atom->shape[i][1] != atom->shape[i][2]))
      error->all("Fix wall/colloid requires spherical particles");

  // insure all particles in group are extended particles

  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int flag = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (atom->shape[type[i]][0] == 0.0) flag = 1;

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall) error->all("Fix wall/colloid requires extended particles");

  FixWall::init();
}

/* ---------------------------------------------------------------------- */

void FixWallColloid::precompute(int m)
{
  coeff1[m] = -4.0/315.0 * epsilon[m] * pow(sigma[m],6.0);
  coeff2[m] = -2.0/3.0 * epsilon[m];
  coeff3[m] = epsilon[m] * pow(sigma[m],6.0)/7560.0;
  coeff4[m] = epsilon[m]/6.0;

  double rinv = 1.0/cutoff[m];
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  offset[m] = coeff3[m]*r4inv*r4inv*rinv - coeff4[m]*r2inv*rinv;
}

/* ----------------------------------------------------------------------
   interaction of all particles in group with all 6 walls (if defined)
   error if any finite-size particle is touching or penetrating wall
------------------------------------------------------------------------- */

void FixWallColloid::wall_particle(int m, double coord)
{
  double delta,delta2,rinv,r2inv,r4inv,r8inv,fwall;
  double r2,rinv2,r2inv2,r4inv2,r6inv2;
  double r3,rinv3,r2inv3,r4inv3,r6inv3;
  double rad,rad2,rad4,rad8,diam,new_coeff2;

  double **x = atom->x;
  double **f = atom->f;
  double **shape = atom->shape;
  int *mask = atom->mask;
  int *type = atom->type;
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
      rad = shape[type[i]][0];
      if (rad >= delta) {
	onflag = 1;
	continue;
      }

      new_coeff2 = coeff2[m]*rad*rad*rad;
      diam = 2.0*rad;
      rad2 = rad*rad;
      rad4 = rad2*rad2;
      rad8 = rad4*rad4;
      delta2 = rad2 - delta*delta;
      rinv = 1.0/delta2;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r8inv = r4inv*r4inv;
      fwall = side * (coeff1[m]*(rad8*rad + 27.0*rad4*rad2*rad*pow(delta,2.0)
				 + 63.0*rad4*rad*pow(delta,4.0)
				 + 21.0*rad2*rad*pow(delta,6.0))*r8inv - 
		      new_coeff2*r2inv);
      f[i][dim] -= fwall;
      r2 = 0.5*diam - delta;
      rinv2 = 1.0/r2;
      r2inv2 = rinv2*rinv2;
      r4inv2 = r2inv2*r2inv2;
      r6inv2 = r4inv2*r2inv2;
      r3 = delta + 0.5*diam;
      rinv3 = 1.0/r3;
      r2inv3 = rinv3*rinv3;
      r4inv3 = r2inv3*r2inv3;
      r6inv3 = r4inv3*r2inv3;
      ewall[0] += coeff3[m]*((-3.5*diam+delta)*r4inv2*r2inv2*rinv2
			     + (3.5*diam+delta)*r4inv3*r2inv3*rinv3) -
	coeff4[m]*((-diam*delta+r2*r3*(log(-r2)-log(r3)))*
		   (-rinv2)*rinv3) - offset[m];
      ewall[m+1] += fwall;
    }

  if (onflag) error->one("Particle on or inside fix wall surface");
}
