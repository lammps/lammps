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
#include "math_extra.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "comm.h"
#include "input.h"
#include "variable.h"
#include "random_mars.h"
#include "memory.h"
#include "error.h"
#include "group.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};
enum{CONSTANT,EQUAL,ATOM};

#define SINERTIA 0.4          // moment of inertia prefactor for sphere
#define EINERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixLangevinEff::FixLangevinEff(LAMMPS *lmp, int narg, char **arg) :
  FixLangevin(lmp, narg, arg)
{
  erforcelangevin = NULL;
}

/* ---------------------------------------------------------------------- */

FixLangevinEff::~FixLangevinEff()
{
  memory->destroy(erforcelangevin);
}

/* ---------------------------------------------------------------------- */

void FixLangevinEff::post_force_no_tally()
{
  double gamma1,gamma2,t_target;

  double **v = atom->v;
  double **f = atom->f;
  double *ervel = atom->ervel;
  double *erforce = atom->erforce;
  int *spin = atom->spin;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double mefactor = domain->dimension/4.0;
  double sqrtmefactor = sqrt(mefactor);

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  // set current t_target and t_sqrt
  // if variable temp, evaluate variable, wrap with clear/add
  // reallocate tforce array if necessary

  if (tstyle == CONSTANT) {
    t_target = t_start + delta * (t_stop-t_start);
    tsqrt = sqrt(t_target);
  } else {
    modify->clearstep_compute();
    if (tstyle == EQUAL) {
      t_target = input->variable->compute_equal(tvar);
      if (t_target < 0.0)
	error->one(FLERR,"Fix langevin/eff variable returned negative temperature");
      tsqrt = sqrt(t_target);
    } else {
      if (nlocal > maxatom2) {
	maxatom2 = atom->nmax;
	memory->destroy(tforce);
	memory->create(tforce,maxatom2,"langevin/eff:tforce");
      }
      input->variable->compute_atom(tvar,igroup,tforce,1,0);
      for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	    if (tforce[i] < 0.0) 
	      error->one(FLERR,
			 "Fix langevin/eff variable returned negative temperature");
    }
    modify->addstep_compute(update->ntimestep + 1);
  }

  // apply damping and thermostat to atoms in group
  // for BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias
  //   test v = 0 since some computes mask non-participating atoms via v = 0
  //   and added force has extra term not multiplied by v = 0
  // for ZEROFLAG:
  //   sum random force over all atoms in group
  //   subtract sum/particles from each atom in group    

  double fran[4],fsum[4],fsumall[4];
  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  int particles = group->count(igroup);
  if (zeroflag) {
    if (particles == 0)
      error->all(FLERR,"Cannot zero Langevin force of 0 atoms/electrons");
  }

  // find number of electrons in group
  int dof,fix_dof;
  dof = domain->dimension * particles;
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);

  // extra_dof = domain->dimension
  dof -= domain->dimension + fix_dof;

  int one = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (fabs(spin[i])==1) one++;
    }
  int nelectrons, dofelectrons, dofnuclei;
  MPI_Allreduce(&one,&nelectrons,1,MPI_INT,MPI_SUM,world);
  dofelectrons = domain->dimension*nelectrons;
  dofnuclei = dof-dofelectrons;

  // thermal partitioning factor between nuclei and electrons 
  // extra dof from electron size
  double gfactor3=(double) (dof+nelectrons)/dofnuclei;

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	gamma1 = gfactor1[type[i]] * gfactor3;
	gamma2 = gfactor2[type[i]] * tsqrt;
	fran[0] = gamma2*(random->uniform()-0.5);
	fran[1] = gamma2*(random->uniform()-0.5);
	fran[2] = gamma2*(random->uniform()-0.5);
	f[i][0] += gamma1*v[i][0] + fran[0];
	f[i][1] += gamma1*v[i][1] + fran[1];
	f[i][2] += gamma1*v[i][2] + fran[2];
	fsum[0] += fran[0];
	fsum[1] += fran[1];
	fsum[2] += fran[2];
        if (fabs(spin[i])==1) {
          fran[3] = sqrtmefactor*gamma2*(random->uniform()-0.5);
          erforce[i] += mefactor*gamma1*ervel[i]+fran[3];
          fsum[3] += fran[3];
        }
      }
    }
  } else if (which == BIAS) {
    temperature->compute_scalar();
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	gamma1 = gfactor1[type[i]] * gfactor3;
	gamma2 = gfactor2[type[i]] * tsqrt;
	temperature->remove_bias(i,v[i]);
	fran[0] = gamma2*(random->uniform()-0.5);
	fran[1] = gamma2*(random->uniform()-0.5);
	fran[2] = gamma2*(random->uniform()-0.5);
	if (v[i][0] != 0.0)
	  f[i][0] += gamma1*v[i][0] + fran[0];
	if (v[i][1] != 0.0)
	  f[i][1] += gamma1*v[i][1] + fran[1];
	if (v[i][2] != 0.0)
	  f[i][2] += gamma1*v[i][2] + fran[2];
	fsum[0] += fran[0];
	fsum[1] += fran[1];
	fsum[2] += fran[2];
        if (fabs(spin[i])==1) {
          fran[3] = sqrtmefactor*gamma2*(random->uniform()-0.5);
          if (ervel[i] != 0.0) erforce[i] += mefactor*gamma1*ervel[i]+fran[3];
          fsum[3] += fran[3];
        }
	temperature->restore_bias(i,v[i]);
      }
    }
  }

  // set total force to zero

  if (zeroflag) {
    MPI_Allreduce(fsum,fsumall,3,MPI_DOUBLE,MPI_SUM,world);
    fsumall[0] /= particles;
    fsumall[1] /= particles;
    fsumall[2] /= particles;
    fsumall[3] /= nelectrons;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	f[i][0] -= fsumall[0];
	f[i][1] -= fsumall[1];
	f[i][2] -= fsumall[2];
        if (fabs(spin[i])==1) erforce[i] -= fsumall[3];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixLangevinEff::post_force_tally()
{
  double gamma1,gamma2,t_target;

  // reallocate flangevin and erforcelangevin if necessary

  if (atom->nlocal > maxatom1) {
    memory->destroy(flangevin);
    memory->destroy(erforcelangevin);
    maxatom1 = atom->nmax;
    memory->create(flangevin,maxatom1,3,"langevin/eff:flangevin");
    memory->create(erforcelangevin,maxatom1,"langevin/eff:erforcelangevin");
  }

  double **v = atom->v;
  double **f = atom->f;
  double *erforce = atom->erforce;
  double *ervel = atom->ervel;
  int *spin = atom->spin;
  double mefactor = domain->dimension/4.0;
  double sqrtmefactor = sqrt(mefactor);

  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  // set current t_target and t_sqrt
  // if variable temp, evaluate variable, wrap with clear/add
  // reallocate tforce array if necessary

  if (tstyle == CONSTANT) {
    t_target = t_start + delta * (t_stop-t_start);
    tsqrt = sqrt(t_target);
  } else {
    modify->clearstep_compute();
    if (tstyle == EQUAL) {
      t_target = input->variable->compute_equal(tvar);
      if (t_target < 0.0)
	error->one(FLERR,"Fix langevin/eff variable returned negative temperature");
      tsqrt = sqrt(t_target);
    } else {
      if (nlocal > maxatom2) {
	maxatom2 = atom->nmax;
	memory->destroy(tforce);
	memory->create(tforce,maxatom2,"langevin/eff:tforce");
      }
      input->variable->compute_atom(tvar,igroup,tforce,1,0);
      for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	    if (tforce[i] < 0.0) 
	      error->one(FLERR,
			 "Fix langevin/eff variable returned negative temperature");
    }
    modify->addstep_compute(update->ntimestep + 1);
  }

  // apply damping and thermostat to appropriate atoms
  // for BIAS:
  //   calculate temperature since some computes require temp
  //   computed on current nlocal atoms to remove bias
  //   test v = 0 since some computes mask non-participating atoms via v = 0
  //   and added force has extra term not multiplied by v = 0

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  int particles = group->count(igroup);
  if (zeroflag) {
    if (particles == 0)
      error->all(FLERR,"Cannot zero Langevin force of 0 atoms/electrons");
  }

  // find number of electrons in group
  int dof,fix_dof;
  dof = domain->dimension * particles;
  fix_dof = 0;
  for (int i = 0; i < modify->nfix; i++)
    fix_dof += modify->fix[i]->dof(igroup);

  // extra_dof = domain->dimension
  dof -= domain->dimension + fix_dof;

  int one = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (fabs(spin[i])==1) one++;
    }
  int nelectrons, dofelectrons, dofnuclei;
  MPI_Allreduce(&one,&nelectrons,1,MPI_INT,MPI_SUM,world);
  dofelectrons = domain->dimension*nelectrons;
  dofnuclei = dof-dofelectrons;

  // thermal partitioning factor between nuclei and electrons
  // extra dof from electron size
  double gfactor3=(double) (dof+nelectrons)/dofnuclei;

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	gamma1 = gfactor1[type[i]] * gfactor3;
	gamma2 = gfactor2[type[i]] * tsqrt;
	flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	f[i][0] += flangevin[i][0];
	f[i][1] += flangevin[i][1];
	f[i][2] += flangevin[i][2];
        if (fabs(spin[i])==1) {
          erforcelangevin[i] = mefactor*gamma1*ervel[i]+sqrtmefactor*gamma2*(random->uniform()-0.5);
          erforce[i] += erforcelangevin[i];
        }
      }
    }
  } else if (which == BIAS) {
    temperature->compute_scalar();
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	gamma1 = gfactor1[type[i]] * gfactor3;
	gamma2 = gfactor2[type[i]] * tsqrt;
	temperature->remove_bias(i,v[i]);
	flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	if (v[i][0] != 0.0) f[i][0] += flangevin[i][0];
	else flangevin[i][0] = 0.0;
	if (v[i][1] != 0.0) f[i][1] += flangevin[i][1];
	else flangevin[i][1] = 0.0;
	if (v[i][2] != 0.0) f[i][2] += flangevin[i][2];
	else flangevin[i][2] = 0.0;
        if (fabs(spin[i])==1) {
          erforcelangevin[i] = mefactor*gamma1*ervel[i]+sqrtmefactor*gamma2*(random->uniform()-0.5);
          if (ervel[i] != 0.0) erforce[i] += erforcelangevin[i];
          else erforcelangevin[i] = 0.0;
        }
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
      if (fabs(spin[i])==1) energy_onestep += erforcelangevin[i];
    }
  energy += energy_onestep*update->dt;
}

/* ---------------------------------------------------------------------- */

double FixLangevinEff::compute_scalar()
{
  if (!tally || flangevin == NULL || erforcelangevin == NULL) return 0.0;

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
        if (fabs(spin[i])==1) energy_onestep += erforcelangevin[i];
      }
    energy = 0.5*energy_onestep*update->dt;
  }

  double energy_me = energy - 0.5*energy_onestep*update->dt;

  double energy_all;	 
  MPI_Allreduce(&energy_me,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);	
  return -energy_all;
}

