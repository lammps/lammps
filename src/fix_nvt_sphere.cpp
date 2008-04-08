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
  if (!atom->omega_flag || !atom->torque_flag)
    error->all("Fix nvt/sphere requires atom attributes omega, torque");

  dttype = new double[atom->ntypes+1];
}

/* ---------------------------------------------------------------------- */

FixNVTSphere::~FixNVTSphere()
{
  delete [] dttype;
}

/* ---------------------------------------------------------------------- */

void FixNVTSphere::init()
{
  FixNVT::init();

  if (!atom->shape)
    error->all("Fix nvt/sphere requires atom attribute shape");

  double **shape = atom->shape;
  for (int i = 1; i <= atom->ntypes; i++)
    if (shape[i][0] != shape[i][1] || shape[i][0] != shape[i][2])
      error->all("Fix nvt/sphere requires spherical particle shapes");
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
  double **quat = atom->quat;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // recompute timesteps since dt may have changed or come via rRESPA

  double dtfrotate = dtf / INERTIA;
  int ntypes = atom->ntypes;
  double **shape = atom->shape;
  for (int i = 1; i <= ntypes; i++)
    dttype[i] = dtfrotate / (0.25*shape[i][0]*shape[i][0]*mass[i]);

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
	
	dtirotate = dttype[itype];
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
	
	dtirotate = dttype[itype];
	omega[i][0] = omega[i][0]*factor + dtirotate*torque[i][0];
	omega[i][1] = omega[i][1]*factor + dtirotate*torque[i][1];
	omega[i][2] = omega[i][2]*factor + dtirotate*torque[i][2];
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
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // recompute timesteps since dt may have changed or come via rRESPA

  double dtfrotate = dtf / INERTIA;
  int ntypes = atom->ntypes;
  double **shape = atom->shape;
  for (int i = 1; i <= ntypes; i++)
    dttype[i] = dtfrotate / (0.25*shape[i][0]*shape[i][0]*mass[i]);

  if (which == NOBIAS) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	itype = type[i];

	dtfm = dtf / mass[itype] * factor;
	v[i][0] = v[i][0]*factor + dtfm*f[i][0];
	v[i][1] = v[i][1]*factor + dtfm*f[i][1];
	v[i][2] = v[i][2]*factor + dtfm*f[i][2];

	dtirotate = dttype[itype];
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
	
	dtirotate = dttype[itype];
	omega[i][0] = (omega[i][0] + dtirotate*torque[i][0]) * factor;
	omega[i][1] = (omega[i][1] + dtirotate*torque[i][1]) * factor;
	omega[i][2] = (omega[i][2] + dtirotate*torque[i][2]) * factor;
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
