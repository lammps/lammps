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

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_nvt_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "comm.h"
#include "group.h"
#include "update.h"
#include "modify.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;

#define INERTIA 0.4          // moment of inertia for sphere

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixNVTSphere::FixNVTSphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVT(lmp, narg, arg)
{
  // error checks

  if (!atom->omega_flag || !atom->torque_flag)
    error->all("Fix nvt/sphere requires atom attributes omega, torque");
  if (!atom->radius_flag && !atom->avec->shape_type)
    error->all("Fix nvt/sphere requires atom attribute radius or shape");
}

/* ---------------------------------------------------------------------- */

void FixNVTSphere::init()
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
	  error->one("Fix nvt/sphere requires extended particles");
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
	  error->one("Fix nvt/sphere requires extended particles");
	if (shape[itype][0] != shape[itype][1] || 
	    shape[itype][0] != shape[itype][2])
	  error->one("Fix nvt/sphere requires spherical particle shapes");
      }
  }

  FixNVT::init();
}

/* ---------------------------------------------------------------------- */

void FixNVTSphere::initial_integrate(int vflag)
{
  int itype;
  double dtfm,dtirotate;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop-t_start);

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
  eta += dtv*eta_dot;
  factor = exp(-dthalf*eta_dot);

  // update v and x of only atoms in group

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
  // 8 cases depending on radius vs shape, rmass vs mass, bias vs nobias

  if (radius) {
    if (rmass) {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    dtfm = dtf / rmass[i];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      } else if (which == BIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / rmass[i];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      }

    } else {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    dtfm = dtf / mass[itype];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*mass[itype]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      } else if (which == BIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / mass[itype];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*mass[itype]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      }
    }

  } else {
    if (rmass) {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    dtfm = dtf / rmass[i];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / (shape[itype][0]*shape[itype][0]*rmass[i]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      } else if (which == BIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / rmass[i];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / (shape[itype][0]*shape[itype][0]*rmass[i]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      }

    } else {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    dtfm = dtf / mass[itype];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / 
	      (shape[itype][0]*shape[itype][0]*mass[itype]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      } else if (which == BIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / mass[itype];
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    x[i][0] += dtv * v[i][0];
	    x[i][1] += dtv * v[i][1];
	    x[i][2] += dtv * v[i][2];
	    
	    dtirotate = dtfrotate / 
	      (shape[itype][0]*shape[itype][0]*mass[itype]);
	    omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	    omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	    omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
	  }
	}
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVTSphere::final_integrate()
{
  int itype;
  double dtfm,dtirotate;

  // update v of only atoms in group

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

  // update v,omega for all particles
  // d_omega/dt = torque / inertia
  // 8 cases depending on radius vs shape, rmass vs mass, bias vs nobias

  if (radius) {
    if (rmass) {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    dtfm = dtf / rmass[i] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      } else {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / rmass[i] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      }

    } else {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    dtfm = dtf / mass[itype] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*mass[itype]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      } else {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / mass[itype] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    
	    dtirotate = dtfrotate / (radius[i]*radius[i]*mass[itype]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      }
    }

  } else {
    if (rmass) {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    dtfm = dtf / rmass[i] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    
	    dtirotate = dtfrotate / (shape[itype][0]*shape[itype][0]*rmass[i]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      } else {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / rmass[i] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    
	    dtirotate = dtfrotate / (shape[itype][0]*shape[itype][0]*rmass[i]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      }

    } else {
      if (which == NOBIAS) {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    dtfm = dtf / mass[itype] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    
	    dtirotate = dtfrotate / 
	      (shape[itype][0]*shape[itype][0]*mass[itype]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      } else {
	for (int i = 0; i < nlocal; i++) {
	  if (mask[i] & groupbit) {
	    itype = type[i];
	    temperature->remove_bias(i,v[i]);
	    dtfm = dtf / mass[itype] * factor;
	    v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	    v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	    v[i][2] = v[i][2]*factor + dtfm*f[i][2];
	    temperature->restore_bias(i,v[i]);
	    
	    dtirotate = dtfrotate / 
	      (shape[itype][0]*shape[itype][0]*mass[itype]);
	    omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	    omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	    omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
	  }
	}
      }
    }
  }

  // compute current T

  t_current = temperature->compute_scalar();

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
}
