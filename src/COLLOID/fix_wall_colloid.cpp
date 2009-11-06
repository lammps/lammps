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
#include "stdlib.h"
#include "string.h"
#include "fix_wall_colloid.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallColloid::FixWallColloid(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 9) error->all("Illegal fix wall/colloid command");

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 3;
  scalar_vector_freq = 1;
  extscalar = 1;
  extvector = 1;

  if (strcmp(arg[3],"xlo") == 0) {
    dim = 0;
    side = -1;
  } else if (strcmp(arg[3],"xhi") == 0) {
    dim = 0;
    side = 1;
  } else if (strcmp(arg[3],"ylo") == 0) {
    dim = 1;
    side = -1;
  } else if (strcmp(arg[3],"yhi") == 0) {
    dim = 1;
    side = 1;
  } else if (strcmp(arg[3],"zlo") == 0) {
    dim = 2;
    side = -1;
  } else if (strcmp(arg[3],"zhi") == 0) {
    dim = 2;
    side = 1;
  } else error->all("Illegal fix wall/colloid command");

  coord = atof(arg[4]);
  epsilon = atof(arg[5]);
  sigma = atof(arg[6]);
  diam = atof(arg[7]);
  cutoff = atof(arg[8]);

  coeff1 = -576.0/315.0 * epsilon * pow(sigma,6.0);
  coeff2 = -288.0/3.0 * 0.125*diam*diam*diam* epsilon;
  coeff3 = 144.0 * epsilon * pow(sigma,6.0)/7560.0;
  coeff4 = 144.0 * epsilon/6.0;

  double rinv = 1.0/cutoff;
  double r2inv = rinv*rinv;
  double r4inv = r2inv*r2inv;
  offset = coeff3*r4inv*r4inv*rinv - coeff4*r2inv*rinv;

  if (dim == 0 && domain->xperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (dim == 1 && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (dim == 2 && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");

  wall_flag = 0;
  wall[0] = wall[1] = wall[2] = wall[3] = 0.0;
}

/* ---------------------------------------------------------------------- */

int FixWallColloid::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallColloid::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallColloid::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallColloid::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallColloid::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta,delta2,rinv,r2inv,r4inv,r8inv,fwall;
  double r2,rinv2,r2inv2,r4inv2,r6inv2;
  double r3,rinv3,r2inv3,r4inv3,r6inv3;
  double rad,rad2,rad4,rad8;
  wall[0] = wall[1] = wall[2] = wall[3] = 0.0;
  wall_flag = 0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side == -1) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta <= 0.0) continue;
      if (delta > cutoff) continue;
      rad = 0.5*diam;
      rad2 = rad*rad;
      rad4 = rad2*rad2;
      rad8 = rad4*rad4;
      delta2 = rad2 - delta*delta;
      rinv = 1.0/delta2;
      r2inv = rinv*rinv;
      r4inv = r2inv*r2inv;
      r8inv = r4inv*r4inv;
      fwall = (coeff1*(rad8*rad + 27.0*rad4*rad2*rad*pow(delta,2.0)
                       + 63.0*rad4*rad*pow(delta,4.0)
                       + 21.0*rad2*rad*pow(delta,6.0))*r8inv
             - coeff2*r2inv) * side;
      f[i][dim] -= fwall;
      r2 = 0.5*diam - delta;
      rinv2 = 1.0/r2;
      r2inv2 = rinv2*rinv2;
      r4inv2 = r2inv2*r2inv2;
      r6inv2 = r4inv2*r2inv2;
      r3 = delta+0.5*diam;
      rinv3 = 1.0/r3;
      r2inv3 = rinv3*rinv3;
      r4inv3 = r2inv3*r2inv3;
      r6inv3 = r4inv3*r2inv3;
      wall[0] += coeff3*((-3.5*diam+delta)*r4inv2*r2inv2*rinv2
			 + (3.5*diam+delta)*r4inv3*r2inv3*rinv3) 
	- coeff4*((-diam*delta+r2*r3*(log(-r2)-log(r3)))*
		  (-rinv2)*rinv3) - offset;
      wall[dim+1] += fwall;
    }
}

/* ---------------------------------------------------------------------- */

void FixWallColloid::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallColloid::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   energy of wall interaction
------------------------------------------------------------------------- */

double FixWallColloid::compute_scalar()
{
  // only sum across procs one time

  if (wall_flag == 0) {
    MPI_Allreduce(wall,wall_all,4,MPI_DOUBLE,MPI_SUM,world);
    wall_flag = 1;
  }
  return wall_all[0];
}

/* ----------------------------------------------------------------------
   components of force on wall
------------------------------------------------------------------------- */

double FixWallColloid::compute_vector(int n)
{
  // only sum across procs one time

  if (wall_flag == 0) {
    MPI_Allreduce(wall,wall_all,4,MPI_DOUBLE,MPI_SUM,world);
    wall_flag = 1;
  }
  return wall_all[n+1];
}
