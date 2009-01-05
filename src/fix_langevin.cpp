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
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "random_mars.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};
enum{MASS,RMASS};

/* ---------------------------------------------------------------------- */

FixLangevin::FixLangevin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal fix langevin command");

  time_depend = 1;

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

  for (int i = 1; i <= atom->ntypes; i++) ratio[i] = 1.0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix langevin command");
      int itype = atoi(arg[iarg+1]);
      double scale = atof(arg[iarg+2]);
      if (itype <= 0 || itype > atom->ntypes)
	error->all("Illegal fix langevin command");
      ratio[itype] = scale;
      iarg += 3;
    } else error->all("Illegal fix langevin command");
  }

  // set temperature = NULL, user can override via fix_modify if wants bias

  id_temp = NULL;
  temperature = NULL;
}

/* ---------------------------------------------------------------------- */

FixLangevin::~FixLangevin()
{
  delete random;
  delete [] gfactor1;
  delete [] gfactor2;
  delete [] ratio;
  delete [] id_temp;
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
  // set force prefactors

  if (atom->mass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      gfactor1[i] = -atom->mass[i] / t_period / force->ftm2v;
      gfactor2[i] = sqrt(atom->mass[i]) * 
	sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) / 
	force->ftm2v;
      gfactor1[i] *= 1.0/ratio[i];
      gfactor2[i] *= 1.0/sqrt(ratio[i]);
    }
  }

  if (temperature && temperature->tempbias) which = BIAS;
  else which = NOBIAS;

  if (atom->mass) massflag = MASS;
  else massflag = RMASS;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixLangevin::setup(int vflag)
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

void FixLangevin::post_force(int vflag)
{
  double gamma1,gamma2;

  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);
  double tsqrt = sqrt(t_target);

  // apply damping and thermostat to appropriate atoms

  if (massflag == MASS) {
    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  gamma1 = gfactor1[type[i]];
	  gamma2 = gfactor2[type[i]] * tsqrt;
	  f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	  f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	  f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	}
      }

    // invoke temperature since some computes require it to remove bias
    // test v = 0 since some computes mask non-participating atoms via v = 0

    } else if (which == BIAS) {
      double tmp = temperature->compute_scalar();
      
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  gamma1 = gfactor1[type[i]];
	  gamma2 = gfactor2[type[i]] * tsqrt;
	  temperature->remove_bias(i,v[i]);
	  if (v[i][0] != 0.0)
	    f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	  if (v[i][1] != 0.0)
	    f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	  if (v[i][2] != 0.0)
	    f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	  temperature->restore_bias(i,v[i]);
	}
      }
    }

  } else {
    double *rmass = atom->rmass;
    double boltz = force->boltz;
    double dt = update->dt;
    double mvv2e = force->mvv2e;
    double ftm2v = force->ftm2v;

    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  gamma1 = -rmass[i] / t_period / ftm2v;
	  gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
	  gamma1 *= 1.0/ratio[type[i]];
	  gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
	  f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	  f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	  f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	}
      }

    // invoke temperature since some computes require it to remove bias
    // test v = 0 since some computes mask non-participating atoms via v = 0

    } else if (which == BIAS) {
      double tmp = temperature->compute_scalar();
      
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  gamma1 = -rmass[i] / t_period / ftm2v;
	  gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
	  gamma1 *= 1.0/ratio[type[i]];
	  gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
	  temperature->remove_bias(i,v[i]);
	  if (v[i][0] != 0.0)
	    f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	  if (v[i][1] != 0.0)
	    f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	  if (v[i][2] != 0.0)
	    f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	  temperature->restore_bias(i,v[i]);
	}
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

/* ---------------------------------------------------------------------- */

void FixLangevin::reset_dt()
{
  if (atom->mass) {
    for (int i = 1; i <= atom->ntypes; i++) {
      gfactor2[i] = sqrt(atom->mass[i]) * 
	sqrt(24.0*force->boltz/t_period/update->dt/force->mvv2e) / 
	force->ftm2v;
      gfactor2[i] *= 1.0/sqrt(ratio[i]);
    }
  }
}

/* ---------------------------------------------------------------------- */

int FixLangevin::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"temp") == 0) {
    if (narg < 2) error->all("Illegal fix_modify command");
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) error->all("Could not find fix_modify temp ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all("Fix_modify temp ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning("Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}
