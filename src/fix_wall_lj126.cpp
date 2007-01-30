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
   Contributing author: Mark Stevens (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_lj126.h"
#include "atom.h"
#include "update.h"
#include "output.h"
#include "respa.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixWallLJ126::FixWallLJ126(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all("Illegal fix wall/lj126 command");

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
  } else error->all("Illegal fix wall/lj126 command");

  coord = atof(arg[4]);
  epsilon = atof(arg[5]);
  sigma = atof(arg[6]);
  cutoff = atof(arg[7]);

  coeff1 = 48.0 * epsilon * pow(sigma,12.0);
  coeff2 = 24.0 * epsilon * pow(sigma,6.0);
  coeff3 = 4.0 * epsilon * pow(sigma,12.0);
  coeff4 = 4.0 * epsilon * pow(sigma,6.0);

  double r2inv = 1.0/(cutoff*cutoff);
  double r6inv = r2inv*r2inv*r2inv;
  offset = r6inv*(coeff3*r6inv - coeff4);
}

/* ---------------------------------------------------------------------- */

int FixWallLJ126::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallLJ126::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallLJ126::setup()
{
  eflag_enable = 1;
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(1);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(1,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
  eflag_enable = 0;
}

/* ---------------------------------------------------------------------- */

void FixWallLJ126::min_setup()
{
  eflag_enable = 1;
  post_force(1);
}

/* ---------------------------------------------------------------------- */

void FixWallLJ126::post_force(int vflag)
{
  bool eflag = false;
  if (eflag_enable) eflag = true;
  else if (output->next_thermo == update->ntimestep) eflag = true;

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta,rinv,r2inv,r6inv,eng;
  if (eflag) eng = 0.0;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (side == -1) delta = x[i][dim] - coord;
      else delta = coord - x[i][dim];
      if (delta <= 0.0) continue;
      if (delta > cutoff) continue;
      rinv = 1.0/delta;
      r2inv = rinv*rinv;
      r6inv = r2inv*r2inv*r2inv;
      f[i][dim] -= r6inv*(coeff1*r6inv - coeff2) * side;
      if (eflag) eng += r6inv*(coeff3*r6inv - coeff4) - offset;
    }

  if (eflag) MPI_Allreduce(&eng,&etotal,1,MPI_DOUBLE,MPI_SUM,world);
}

/* ---------------------------------------------------------------------- */

void FixWallLJ126::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixWallLJ126::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

double FixWallLJ126::thermo(int n)
{
  if (n == 0) return etotal;
  else return 0.0;
}
