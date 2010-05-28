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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "fix_nh_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "group.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia for sphere

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixNHSphere::FixNHSphere(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  if (!atom->omega_flag || !atom->torque_flag)
    error->all("Fix nvt/nph/npt sphere requires "
	       "atom attributes omega, torque");
  if (!atom->radius_flag && !atom->avec->shape_type)
    error->all("Fix nvt/nph/npt sphere requires "
	       "atom attribute radius or shape");
}

/* ---------------------------------------------------------------------- */

void FixNHSphere::init()
{
  int i,itype;

  // check that all particles are finite-size and spherical
  // no point particles allowed

  if (atom->radius_flag) {
    double *radius = atom->radius;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	if (radius[i] == 0.0)
	  error->one("Fix nvt/sphere requires extended particles");
      }

  } else {
    double **shape = atom->shape;
    int *type = atom->type;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;

    for (i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
	itype = type[i];
	if (shape[itype][0] == 0.0)
	  error->one("Fix nvt/sphere requires extended particles");
	if (shape[itype][0] != shape[itype][1] || 
	    shape[itype][0] != shape[itype][2])
	  error->one("Fix nvt/sphere requires spherical particle shapes");
      }
  }

  FixNH::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of velocities 
-----------------------------------------------------------------------*/

void FixNHSphere::nve_v()
{
  // standard nve_v velocity update

  FixNH::nve_v();

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
  double dtirotate;
  int itype;

  // update omega for all particles
  // d_omega/dt = torque / inertia
  // 4 cases depending on radius vs shape and rmass vs mass

  if (radius) {
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {    
	if (mask[i] & groupbit) {
	  dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
	  omega[i][0] += dtirotate*torque[i][0];
	  omega[i][1] += dtirotate*torque[i][1];
	  omega[i][2] += dtirotate*torque[i][2];
	}
      }
    } else {
      for (int i = 0; i < nlocal; i++) {    
	if (mask[i] & groupbit) {
	  dtirotate = dtfrotate / (radius[i]*radius[i]*mass[type[i]]);
	  omega[i][0] += dtirotate*torque[i][0];
	  omega[i][1] += dtirotate*torque[i][1];
	  omega[i][2] += dtirotate*torque[i][2];
	}
      }
    }

  } else {
    if (rmass) {
      for (int i = 0; i < nlocal; i++) {    
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtirotate = dtfrotate / (shape[itype][0]*shape[itype][0]*rmass[i]);
	  omega[i][0] += dtirotate*torque[i][0];
	  omega[i][1] += dtirotate*torque[i][1];
	  omega[i][2] += dtirotate*torque[i][2];
	}
      }
    } else {
      for (int i = 0; i < nlocal; i++) {    
	if (mask[i] & groupbit) {
	  itype = type[i];
	  dtirotate = dtfrotate / 
	    (shape[itype][0]*shape[itype][0]*mass[itype]);
	  omega[i][0] += dtirotate*torque[i][0];
	  omega[i][1] += dtirotate*torque[i][1];
	  omega[i][2] += dtirotate*torque[i][2];
	}
      }
    }
  }

}

/* ----------------------------------------------------------------------
   perform half-step scaling of velocities
-----------------------------------------------------------------------*/

void FixNHSphere::nh_v_temp()
{
  // standard nh_v_temp velocity update

  FixNH::nh_v_temp();

  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double factor_rotate = exp(-dthalf*eta_dot[0]);

  for (int i = 0; i < nlocal; i++) {    
    if (mask[i] & groupbit) {
      omega[i][0] *= factor_rotate;
      omega[i][1] *= factor_rotate;
      omega[i][2] *= factor_rotate;
    }
  }
}
