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
   Contributing author: Andres Jaramillo-Botero
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_langevin_eff.h"
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
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixLangevinEff::FixLangevinEff(LAMMPS *lmp, int narg, char **arg) :
  FixLangevin(lmp, narg, arg)
{
  erforcelangevin = NULL;
}

/* ---------------------------------------------------------------------- */

FixLangevinEff::~FixLangevinEff()
{
  memory->sfree(erforcelangevin);
}

/* ---------------------------------------------------------------------- */

void FixLangevinEff::post_force_no_tally()
{
  double gamma1,gamma2;

  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *ervel = atom->ervel;
  double *erforce = atom->erforce;
  int *spin = atom->spin;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);
  double tsqrt = sqrt(t_target);

  // apply damping and thermostat to atoms in group
  // for BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias
  //   test v = 0 since some computes mask non-participating atoms via v = 0
  //   and added force has extra term not multiplied by v = 0

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        gamma1 = gfactor1[type[i]];
	gamma2 = gfactor2[type[i]] * tsqrt;
	f[i][0] += gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	f[i][1] += gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	f[i][2] += gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
        if (abs(spin[i])==1) erforce[i] += 0.75*gamma1*ervel[i] + 0.866025404*gamma2*(random->uniform()-0.5);
      }
    }
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
        if (abs(spin[i])==1 && ervel[i] != 0.0)
          erforce[i] += 0.75*gamma1*ervel[i] + 0.866025404*gamma2*(random->uniform()-0.5);
	temperature->restore_bias(i,v[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLangevinEff::post_force_tally()
{
  double gamma1,gamma2;

  // reallocate flangevin if necessary

  if (atom->nmax > nmax) {
    memory->destroy(flangevin);
    memory->sfree(erforcelangevin);
    nmax = atom->nmax;
    memory->create(flangevin,nmax,3,"langevin:flangevin");
    erforcelangevin = (double *) 
      memory->smalloc(nmax*sizeof(double),"langevin/eff:erforcelangevin");
  }

  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *erforce = atom->erforce;
  double *ervel = atom->ervel;
  int *spin = atom->spin;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  double t_target = t_start + delta * (t_stop-t_start);
  double tsqrt = sqrt(t_target);

  // apply damping and thermostat to appropriate atoms
  // for BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias
  //   test v = 0 since some computes mask non-participating atoms via v = 0
  //   and added force has extra term not multiplied by v = 0

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	gamma1 = gfactor1[type[i]];
	gamma2 = gfactor2[type[i]] * tsqrt;
	flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
        erforcelangevin[i] = 0.75*gamma1*ervel[i]+0.866025404*gamma2*(random->uniform()-0.5);
	f[i][0] += flangevin[i][0];
	f[i][1] += flangevin[i][1];
	f[i][2] += flangevin[i][2];
        if (abs(spin[i])==1) erforce[i] += erforcelangevin[i];
      }
    }
  } else if (which == BIAS) {
    double tmp = temperature->compute_scalar();
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	gamma1 = gfactor1[type[i]];
	gamma2 = gfactor2[type[i]] * tsqrt;
	temperature->remove_bias(i,v[i]);
	flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
        erforcelangevin[i] = 0.75*gamma1*ervel[i]+0.866025404*gamma2*(random->uniform()-0.5);
	if (v[i][0] != 0.0) f[i][0] += flangevin[i][0];
	else flangevin[i][0] = 0.0;
	if (v[i][1] != 0.0) f[i][1] += flangevin[i][1];
	else flangevin[i][1] = 0.0;
	if (v[i][2] != 0.0) f[i][2] += flangevin[i][2];
	else flangevin[i][2] = 0.0;
        if (abs(spin[i])==1 && ervel[i] != 0.0) erforce[i] += erforcelangevin[i];
	temperature->restore_bias(i,v[i]);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy transfer to thermal reservoir
------------------------------------------------------------------------- */

void FixLangevinEff::end_of_step()
{
  if (!tally) return;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;  
  int *spin = atom->spin;

  energy_onestep = 0.0;
 
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      energy_onestep += flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] + 
	flangevin[i][2]*v[i][2];
      if (abs(spin[i])==1) energy_onestep += erforcelangevin[i];
    }
  energy += energy_onestep*update->dt;
}

/* ---------------------------------------------------------------------- */

double FixLangevinEff::compute_scalar()
{
  if (!tally) return 0.0;

  // capture the very first energy transfer to thermal reservoir

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;  
  int *spin = atom->spin;

  if (update->ntimestep == update->beginstep) {
    energy_onestep = 0.0;
    for (int i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) {
	energy_onestep += flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] + 
	  flangevin[i][2]*v[i][2];
        if (abs(spin[i])==1) energy_onestep += erforcelangevin[i];
      }
    energy = 0.5*energy_onestep*update->dt;
  }

  double energy_me = energy - 0.5*energy_onestep*update->dt;

  double energy_all;	 
  MPI_Allreduce(&energy_me,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);	
  return -energy_all;
}
