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

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_langevin.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixLangevin::FixLangevin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal fix langevin command");

  t_start = atof(arg[3]);
  t_stop = atof(arg[4]);
  t_period = atof(arg[5]);
  int seed = atoi(arg[6]);

  if (t_period <= 0.0) error->all("Fix langevin period must be > 0.0");
  if (seed <= 0) error->all("Illegal fix langevin command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];
  ratio = new double[atom->ntypes+1];
  
  // optional args

  flagx = flagy = flagz = 1;
  for (int i = 1; i <= atom->ntypes; i++) ratio[i] = 1.0;
  iregion = -1;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"axes") == 0) {
      if (iarg+4 > narg) error->all("Illegal fix langevin command");
      flagx = atoi(arg[iarg+1]);
      flagy = atoi(arg[iarg+2]);
      flagz = atoi(arg[iarg+3]);
      iarg += 4;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix langevin command");
      int itype = atoi(arg[iarg+1]);
      double scale = atof(arg[iarg+2]);
      if (itype <= 0 || itype > atom->ntypes)
	error->all("Illegal fix langevin command");
      ratio[itype] = scale;
      iarg += 3;
    } else if (strcmp(arg[iarg],"region") == 0) {
      iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all("Fix langevin region ID does not exist");
      iarg += 2;
    } else error->all("Illegal fix langevin command");
  }
}

/* ---------------------------------------------------------------------- */

FixLangevin::~FixLangevin()
{
  delete random;
  delete [] gfactor1;
  delete [] gfactor2;
  delete [] ratio;
}

/* ---------------------------------------------------------------------- */

int FixLangevin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevin::init()
{
  if (atom->mass == NULL)
    error->all("Cannot use fix langevin without per-type mass defined");

  // set force prefactors

  for (int i = 1; i <= atom->ntypes; i++) {
    gfactor1[i] = - atom->mass[i] / t_period / force->ftm2v;
    gfactor2[i] = sqrt(atom->mass[i]) * 
      sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) / force->ftm2v;
    gfactor1[i] *= 1.0/ratio[i];
    gfactor2[i] *= 1.0/sqrt(ratio[i]);
  }

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixLangevin::setup()
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

void FixLangevin::post_force(int vflag)
{
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);
  double tsqrt = sqrt(t_target);

  double gamma1,gamma2;

  // apply damping and thermostat to all atoms in fix group

  if (iregion == -1) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	gamma1 = gfactor1[type[i]];
	gamma2 = gfactor2[type[i]] * tsqrt;
	if (flagx) f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	if (flagy) f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	if (flagz) f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
      }
    }

  // apply damping and thermostat to all atoms in fix group and in region

  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit &&
	  domain->regions[iregion]->match(x[i][0],x[i][1],x[i][2])) {
	gamma1 = gfactor1[type[i]];
	gamma2 = gfactor2[type[i]] * tsqrt;
	if (flagx) f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	if (flagy) f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	if (flagz) f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLangevin::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLangevin::reset_target(double t_new)
{
  t_start = t_stop = t_new;
}
