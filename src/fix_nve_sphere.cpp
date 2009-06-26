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
#include "stdio.h"
#include "string.h"
#include "fix_nve_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia for sphere

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

FixNVESphere::FixNVESphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg < 3) error->all("Illegal fix nve/sphere command");

  time_integrate = 1;

  // process extra keywords

  extra = NONE;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"update") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix nve/sphere command");
      if (strcmp(arg[iarg+1],"dipole") == 0) extra = DIPOLE;
      else error->all("Illegal fix nve/sphere command");
      iarg += 2;
    } else error->all("Illegal fix nve/sphere command");
  }

  // error checks

  if (!atom->omega_flag || !atom->torque_flag)
    error->all("Fix nve/sphere requires atom attributes omega, torque");
  if (!atom->radius_flag && !atom->avec->shape_type)
    error->all("Fix nve/sphere requires atom attribute radius or shape");
  if (extra == DIPOLE && !atom->mu_flag)
    error->all("Fix nve/sphere requires atom attribute mu");
}

/* ---------------------------------------------------------------------- */

int FixNVESphere::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::init()
{
  int i,itype;

  // check that all particles are finite-size and spherical
  // no point particles allowed

  if (atom->radius_flag) {
    double *radius = atom->radius;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (radius[i] == 0.0)
	  error->one("Fix nve/sphere requires extended particles");
      }

  } else {
    double **shape = atom->shape;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    if (igroup == atom->firstgroup) nlocal = atom->nfirst;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	itype = type[i];
	if (shape[itype][0] == 0.0)
	  error->one("Fix nve/sphere requires extended particles");
	if (shape[itype][0] != shape[itype][1] || 
	    shape[itype][0] != shape[itype][2])
	  error->one("Fix nve/sphere requires spherical particle shapes");
      }
  }

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::initial_integrate(int vflag)
{
  int itype;
  double dtfm,dtirotate,msq,scale;
  double g[3];

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / INERTIA;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia
  // 4 cases depending on radius vs shape and rmass vs mass

  if (radius) {
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  x[i][0] += dtv * v[i][0];
	  x[i][1] += dtv * v[i][1];
	  x[i][2] += dtv * v[i][2];
	  
	  dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }

    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtfm = dtf / mass[itype];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  x[i][0] += dtv * v[i][0];
	  x[i][1] += dtv * v[i][1];
	  x[i][2] += dtv * v[i][2];
	  
	  dtirotate = dtfrotate / (radius[i]*radius[i]*mass[itype]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }
    }

  } else {
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  x[i][0] += dtv * v[i][0];
	  x[i][1] += dtv * v[i][1];
	  x[i][2] += dtv * v[i][2];
	  
	  dtirotate = dtfrotate / (shape[itype][0]*shape[itype][0]*rmass[i]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }

    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtfm = dtf / mass[itype];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  x[i][0] += dtv * v[i][0];
	  x[i][1] += dtv * v[i][1];
	  x[i][2] += dtv * v[i][2];
	  
	  dtirotate = dtfrotate / 
	    (shape[itype][0]*shape[itype][0]*mass[itype]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }
    }
  }

  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  if (extra == DIPOLE) {
    double **mu = atom->mu;
    double *dipole = atom->dipole;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	if (dipole[type[i]] > 0.0) {
	  g[0] = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
	  g[1] = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
	  g[2] = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
	  msq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
	  scale = dipole[type[i]]/sqrt(msq);
	  mu[i][0] = g[0]*scale;
	  mu[i][1] = g[1]*scale;
	  mu[i][2] = g[2]*scale;
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphere::final_integrate()
{
  int itype;
  double dtfm,dtirotate;

  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  double *radius = atom->radius;
  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / INERTIA;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia
  // 4 cases depending on radius vs shape and rmass vs mass

  if (radius) {
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  
	  dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }

    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtfm = dtf / mass[itype];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  
	  dtirotate = dtfrotate / (radius[i]*radius[i]*mass[itype]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }
    }

  } else {
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtfm = dtf / rmass[i];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  
	  dtirotate = dtfrotate / (shape[itype][0]*shape[itype][0]*rmass[i]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }

    } else {
      for (int i = 0; i < nlocal; i++) {
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtfm = dtf / mass[itype];
	  v[i][0] += dtfm * f[i][0];
	  v[i][1] += dtfm * f[i][1];
	  v[i][2] += dtfm * f[i][2];
	  
	  dtirotate = dtfrotate / 
	    (shape[itype][0]*shape[itype][0]*mass[itype]);
	  omega[i][0] += dtirotate * torque[i][0];
	  omega[i][1] += dtirotate * torque[i][1];
	  omega[i][2] += dtirotate * torque[i][2];
	}
      }
    }
  }
}
