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
   Contributing author: Carolyn Phillips (U Mich), reservoir energy tally
------------------------------------------------------------------------- */

#include "mpi.h"
#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_langevin.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
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

FixLangevin::FixLangevin(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all(FLERR,"Illegal fix langevin command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  nevery = 1;

  tstr = NULL;
  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    tstr = new char[n];
    strcpy(tstr,&arg[3][2]);
  } else {
    t_start = atof(arg[3]);
    tstyle = CONSTANT;
  }

  t_stop = atof(arg[4]);
  t_period = atof(arg[5]);
  int seed = atoi(arg[6]);

  if (t_period <= 0.0) error->all(FLERR,"Fix langevin period must be > 0.0");
  if (seed <= 0) error->all(FLERR,"Illegal fix langevin command");

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

  // allocate per-type arrays for force prefactors

  gfactor1 = new double[atom->ntypes+1];
  gfactor2 = new double[atom->ntypes+1];
  ratio = new double[atom->ntypes+1];
  
  // optional args

  for (int i = 1; i <= atom->ntypes; i++) ratio[i] = 1.0;
  oflag = aflag = 0;
  tally = 0;
  zeroflag = 0;

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"angmom") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix langevin command");
      if (strcmp(arg[iarg+1],"no") == 0) aflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) aflag = 1;
      else error->all(FLERR,"Illegal fix langevin command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix langevin command");
      if (strcmp(arg[iarg+1],"no") == 0) oflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) oflag = 1;
      else error->all(FLERR,"Illegal fix langevin command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"scale") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix langevin command");
      int itype = atoi(arg[iarg+1]);
      double scale = atof(arg[iarg+2]);
      if (itype <= 0 || itype > atom->ntypes)
	error->all(FLERR,"Illegal fix langevin command");
      ratio[itype] = scale;
      iarg += 3;
    } else if (strcmp(arg[iarg],"tally") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix langevin command");
      if (strcmp(arg[iarg+1],"no") == 0) tally = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) tally = 1;
      else error->all(FLERR,"Illegal fix langevin command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"zero") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix langevin command");
      if (strcmp(arg[iarg+1],"no") == 0) zeroflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) zeroflag = 1;
      else error->all(FLERR,"Illegal fix langevin command");
      iarg += 2;
    } else error->all(FLERR,"Illegal fix langevin command");
  }

  // error check

  if (aflag) {
    avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
    if (!avec) 
      error->all(FLERR,"Fix langevin angmom requires atom style ellipsoid");
  }

  // set temperature = NULL, user can override via fix_modify if wants bias

  id_temp = NULL;
  temperature = NULL;

  // flangevin is unallocated until first call to setup()
  // compute_scalar checks for this and returns 0.0 if flangevin is NULL

  energy = 0.0;
  flangevin = NULL;
  tforce = NULL;
  maxatom1 = maxatom2 = 0;
}

/* ---------------------------------------------------------------------- */

FixLangevin::~FixLangevin()
{
  delete random;
  delete [] tstr;
  delete [] gfactor1;
  delete [] gfactor2;
  delete [] ratio;
  delete [] id_temp;
  memory->destroy(flangevin);
  memory->destroy(tforce);
}

/* ---------------------------------------------------------------------- */

int FixLangevin::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= END_OF_STEP;
  mask |= THERMO_ENERGY;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLangevin::init()
{
  if (oflag && !atom->sphere_flag)
    error->all(FLERR,"Fix langevin omega requires atom style sphere");
  if (aflag && !atom->ellipsoid_flag)
    error->all(FLERR,"Fix langevin angmom requires atom style ellipsoid");

  // check variable

  if (tstr) {
    tvar = input->variable->find(tstr);
    if (tvar < 0) 
      error->all(FLERR,"Variable name for fix langevin does not exist");
    if (input->variable->equalstyle(tvar)) tstyle = EQUAL;
    else if (input->variable->atomstyle(tvar)) tstyle = ATOM;
    else error->all(FLERR,"Variable for fix langevin is invalid style");
  }

  // if oflag or aflag set, check that all group particles are finite-size

  if (oflag) {
    double *radius = atom->radius;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	if (radius[i] == 0.0)
	  error->one(FLERR,"Fix langevin omega requires extended particles");
  }

  if (aflag) {
    int *ellipsoid = atom->ellipsoid;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
	if (ellipsoid[i] < 0)
	  error->one(FLERR,"Fix langevin angmom requires extended particles");
  }

  // set force prefactors

  if (!atom->rmass) {
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

  if (strstr(update->integrate_style,"respa"))
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixLangevin::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
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
  if (tally) post_force_tally();
  else post_force_no_tally();
}

/* ---------------------------------------------------------------------- */

void FixLangevin::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixLangevin::post_force_no_tally()
{
  double gamma1,gamma2,t_target;

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
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
	error->one(FLERR,"Fix langevin variable returned negative temperature");
      tsqrt = sqrt(t_target);
    } else {
      if (nlocal > maxatom2) {
	maxatom2 = atom->nmax;
	memory->destroy(tforce);
	memory->create(tforce,maxatom2,"langevin:tforce");
      }
      input->variable->compute_atom(tvar,igroup,tforce,1,0);
      for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	    if (tforce[i] < 0.0) 
	      error->one(FLERR,
			 "Fix langevin variable returned negative temperature");
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
  //   subtract sum/count from each atom in group    

  double fran[3],fsum[3],fsumall[3];
  fsum[0] = fsum[1] = fsum[2] = 0.0;
  bigint count;

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;
  
  if (zeroflag) {
    count = group->count(igroup);
    if (count == 0)
      error->all(FLERR,"Cannot zero Langevin force of 0 atoms");
  }
  
  if (rmass) {
    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = -rmass[i] / t_period / ftm2v;
	  gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
	  gamma1 *= 1.0/ratio[type[i]];
	  gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
	  fran[0] = gamma2*(random->uniform()-0.5);
	  fran[1] = gamma2*(random->uniform()-0.5);
	  fran[2] = gamma2*(random->uniform()-0.5);
	  f[i][0] += gamma1*v[i][0] + fran[0];
	  f[i][1] += gamma1*v[i][1] + fran[1];
	  f[i][2] += gamma1*v[i][2] + fran[2];
	  fsum[0] += fran[0];
	  fsum[1] += fran[1];
	  fsum[2] += fran[2];
	}
      }

    } else if (which == BIAS) {
      temperature->compute_scalar();
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = -rmass[i] / t_period / ftm2v;
	  gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
	  gamma1 *= 1.0/ratio[type[i]];
	  gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
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
	  temperature->restore_bias(i,v[i]);
	}
      }
    }

  } else {
    
    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = gfactor1[type[i]];
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
	}
      }

    } else if (which == BIAS) {
      temperature->compute_scalar();
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = gfactor1[type[i]];
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
	  temperature->restore_bias(i,v[i]);
	}
      }
    }
  }

  // set total force to zero

  if (zeroflag) {
    MPI_Allreduce(fsum,fsumall,3,MPI_DOUBLE,MPI_SUM,world);
    fsumall[0] /= count;
    fsumall[1] /= count;
    fsumall[2] /= count;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	f[i][0] -= fsumall[0];
	f[i][1] -= fsumall[1];
	f[i][2] -= fsumall[2];
      }
    }
  }

  // thermostat omega and angmom

  if (oflag) omega_thermostat();
  if (aflag) angmom_thermostat();
}

/* ---------------------------------------------------------------------- */

void FixLangevin::post_force_tally()
{
  double gamma1,gamma2,t_target;

  // reallocate flangevin if necessary

  if (atom->nlocal > maxatom1) {
    memory->destroy(flangevin);
    maxatom1 = atom->nmax;
    memory->create(flangevin,maxatom1,3,"langevin:flangevin");
  }

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
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
	error->one(FLERR,"Fix langevin variable returned negative temperature");
      tsqrt = sqrt(t_target);
    } else {
      if (nlocal > maxatom2) {
	maxatom2 = atom->nmax;
	memory->destroy(tforce);
	memory->create(tforce,maxatom2,"langevin:tforce");
      }
      input->variable->compute_atom(tvar,igroup,tforce,1,0);
      for (int i = 0; i < nlocal; i++)
	if (mask[i] & groupbit)
	    if (tforce[i] < 0.0) 
	      error->one(FLERR,
			 "Fix langevin variable returned negative temperature");
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

  if (rmass) {
    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = -rmass[i] / t_period / ftm2v;
	  gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
	  gamma1 *= 1.0/ratio[type[i]];
	  gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
	  flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	  flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	  flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	  f[i][0] += flangevin[i][0];
	  f[i][1] += flangevin[i][1];
	  f[i][2] += flangevin[i][2];
	}
      }

    } else if (which == BIAS) {
      temperature->compute_scalar();
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = -rmass[i] / t_period / ftm2v;
	  gamma2 = sqrt(rmass[i]) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
	  gamma1 *= 1.0/ratio[type[i]];
	  gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
	  temperature->remove_bias(i,v[i]);
	  flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	  flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	  flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	  if (v[i][0] != 0.0) f[i][0] += flangevin[i][0];
	  else flangevin[i][0] = 0;
	  if (v[i][1] != 0.0) f[i][1] += flangevin[i][1];
	  else flangevin[i][1] = 0;
	  if (v[i][2] != 0.0) f[i][2] += flangevin[i][2];
	  else flangevin[i][2] = 0;
	  temperature->restore_bias(i,v[i]);
	}
      }
    }

  } else {
    if (which == NOBIAS) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = gfactor1[type[i]];
	  gamma2 = gfactor2[type[i]] * tsqrt;
	  flangevin[i][0] = gamma1*v[i][0] + gamma2*(random->uniform()-0.5);
	  flangevin[i][1] = gamma1*v[i][1] + gamma2*(random->uniform()-0.5);
	  flangevin[i][2] = gamma1*v[i][2] + gamma2*(random->uniform()-0.5);
	  f[i][0] += flangevin[i][0];
	  f[i][1] += flangevin[i][1];
	  f[i][2] += flangevin[i][2];
	}
      }

    } else if (which == BIAS) {
      temperature->compute_scalar();
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
	  gamma1 = gfactor1[type[i]];
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
	  temperature->restore_bias(i,v[i]);
	}
      }
    }
  }

  // thermostat omega and angmom

  if (oflag) omega_thermostat();
  if (aflag) angmom_thermostat();
}

/* ----------------------------------------------------------------------
   thermostat rotational dof via omega
------------------------------------------------------------------------- */

void FixLangevin::omega_thermostat()
{
  double gamma1,gamma2;

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  double **torque = atom->torque;
  double **omega = atom->omega;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double tran[3];
  double inertiaone;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      inertiaone = SINERTIA*radius[i]*radius[i]*rmass[i];
      if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
      gamma1 = -inertiaone / t_period / ftm2v;
      gamma2 = sqrt(inertiaone) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
      gamma1 *= 1.0/ratio[type[i]];
      gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
      tran[0] = gamma2*(random->uniform()-0.5);
      tran[1] = gamma2*(random->uniform()-0.5);
      tran[2] = gamma2*(random->uniform()-0.5);
      torque[i][0] += gamma1*omega[i][0] + tran[0];
      torque[i][1] += gamma1*omega[i][1] + tran[1];
      torque[i][2] += gamma1*omega[i][2] + tran[2];
    }
  }
}

/* ----------------------------------------------------------------------
   thermostat rotational dof via angmom
------------------------------------------------------------------------- */

void FixLangevin::angmom_thermostat()
{
  double gamma1,gamma2;

  double boltz = force->boltz;
  double dt = update->dt;
  double mvv2e = force->mvv2e;
  double ftm2v = force->ftm2v;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  double **torque = atom->torque;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  double inertia[3],omega[3],tran[3];
  double *shape,*quat;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      shape = bonus[ellipsoid[i]].shape;
      inertia[0] = EINERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      inertia[1] = EINERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      inertia[2] = EINERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
      quat = bonus[ellipsoid[i]].quat;
      MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
      
      if (tstyle == ATOM) tsqrt = sqrt(tforce[i]);
      gamma1 = -1.0 / t_period / ftm2v;
      gamma2 = sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v;
      gamma1 *= 1.0/ratio[type[i]];
      gamma2 *= 1.0/sqrt(ratio[type[i]]) * tsqrt;
      tran[0] = sqrt(inertia[0])*gamma2*(random->uniform()-0.5);
      tran[1] = sqrt(inertia[1])*gamma2*(random->uniform()-0.5);
      tran[2] = sqrt(inertia[2])*gamma2*(random->uniform()-0.5);
      torque[i][0] += inertia[0]*gamma1*omega[0] + tran[0];
      torque[i][1] += inertia[1]*gamma1*omega[1] + tran[1];
      torque[i][2] += inertia[2]*gamma1*omega[2] + tran[2];
    }
  }
}

/* ----------------------------------------------------------------------
   tally energy transfer to thermal reservoir
------------------------------------------------------------------------- */

void FixLangevin::end_of_step()
{
  if (!tally) return;

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;  

  energy_onestep = 0.0;
 
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      energy_onestep += flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] + 
	flangevin[i][2]*v[i][2];

  energy += energy_onestep*update->dt;
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
    if (narg < 2) error->all(FLERR,"Illegal fix_modify command");
    delete [] id_temp;
    int n = strlen(arg[1]) + 1;
    id_temp = new char[n];
    strcpy(id_temp,arg[1]);

    int icompute = modify->find_compute(id_temp);
    if (icompute < 0) 
      error->all(FLERR,"Could not find fix_modify temperature ID");
    temperature = modify->compute[icompute];

    if (temperature->tempflag == 0)
      error->all(FLERR,
		 "Fix_modify temperature ID does not compute temperature");
    if (temperature->igroup != igroup && comm->me == 0)
      error->warning(FLERR,"Group for fix_modify temp != fix group");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

double FixLangevin::compute_scalar()
{
  if (!tally || flangevin == NULL) return 0.0;

  // capture the very first energy transfer to thermal reservoir

  double **v = atom->v;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;  

  if (update->ntimestep == update->beginstep) {
    energy_onestep = 0.0;
    for (int i = 0; i < nlocal; i++) 
      if (mask[i] & groupbit) 
	energy_onestep += flangevin[i][0]*v[i][0] + flangevin[i][1]*v[i][1] + 
	  flangevin[i][2]*v[i][2];
    energy = 0.5*energy_onestep*update->dt;
  }

  double energy_me = energy - 0.5*energy_onestep*update->dt;

  double energy_all;	 
  MPI_Allreduce(&energy_me,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);	
  return -energy_all;
}

/* ----------------------------------------------------------------------
   memory usage of tally array
------------------------------------------------------------------------- */

double FixLangevin::memory_usage()
{
  double bytes = 0.0;
  if (tally) double bytes = atom->nmax*3 * sizeof(double);
  if (tforce) bytes = atom->nmax * sizeof(double);
  return bytes;
}
