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

#include "string.h"
#include "stdlib.h"
#include "math.h"
#include "fix_npt_asphere.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "compute.h"
#include "kspace.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixNPTAsphere::FixNPTAsphere(LAMMPS *lmp, int narg, char **arg) :
  FixNPT(lmp, narg, arg)
{
  inertia = 
    memory->create_2d_double_array(atom->ntypes+1,3,"fix_npt_asphere:inertia");

  // error checks

  if (!atom->quat_flag || !atom->angmom_flag || !atom->torque_flag ||
      !atom->avec->shape_type)
    error->all("Fix npt/asphere requires atom attributes "
	       "quat, angmom, torque, shape");
  if (atom->radius_flag || atom->rmass_flag)
    error->all("Fix npt/asphere cannot be used with atom attributes "
	       "diameter or rmass");
}

/* ---------------------------------------------------------------------- */

FixNPTAsphere::~FixNPTAsphere()
{
  memory->destroy_2d_double_array(inertia);
}

/* ---------------------------------------------------------------------- */

void FixNPTAsphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (shape[type[i]][0] == 0.0)
	error->one("Fix nvt/asphere requires extended particles");

  FixNPT::init();
  calculate_inertia();
}

/* ---------------------------------------------------------------------- */

void FixNPTAsphere::initial_integrate(int vflag)
{
  int i;
  double dtfm;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;

  // update eta_dot

  t_target = t_start + delta * (t_stop-t_start);
  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
  eta += dtv*eta_dot;

  // update omega_dot
  // for non-varying dims, p_freq is 0.0, so omega_dot doesn't change

  double f_omega,volume;
  if (dimension == 3) volume = domain->xprd*domain->yprd*domain->zprd;
  else volume = domain->xprd*domain->yprd;
  double denskt = atom->natoms*boltz*t_target / volume * nktv2p;

  for (i = 0; i < 3; i++) {
    p_target[i] = p_start[i] + delta * (p_stop[i]-p_start[i]);
    f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
    omega_dot[i] += f_omega*dthalf;
    omega_dot[i] *= drag_factor;
    omega[i] += dtv*omega_dot[i];
    factor[i] = exp(-dthalf*(eta_dot+omega_dot[i]));
    dilation[i] = exp(dthalf*omega_dot[i]);
  }
  factor_rotate = exp(-dthalf*eta_dot);

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **quat = atom->quat;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (which == NOBIAS) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
	v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
	v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
      }
    }
  } else if (which == BIAS) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	temperature->remove_bias(i,v[i]);
	dtfm = dtf / mass[type[i]];
	v[i][0] = v[i][0]*factor[0] + dtfm*f[i][0];
	v[i][1] = v[i][1]*factor[1] + dtfm*f[i][1];
	v[i][2] = v[i][2]*factor[2] + dtfm*f[i][2];
	temperature->restore_bias(i,v[i]);
      }
    }
  }

  // remap simulation box and all owned atoms by 1/2 step

  remap(0);

  // x update by full step only for atoms in group

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];
    }
  }

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  // update angular momentum by 1/2 step for all particles
  // update quaternion a full step via Richardson iteration
  // returns new normalized quaternion

  for (i = 0; i < nlocal; i++) {    
    if (mask[i] & groupbit) {
      angmom[i][0] = angmom[i][0]*factor_rotate + dtf*torque[i][0];
      angmom[i][1] = angmom[i][1]*factor_rotate + dtf*torque[i][1];
      angmom[i][2] = angmom[i][2]*factor_rotate + dtf*torque[i][2];
		
      richardson(quat[i],angmom[i],inertia[type[i]]);
    }
  }

  // remap simulation box and all owned atoms by 1/2 step
  // redo KSpace coeffs since volume has changed

  remap(0);
  if (kspace_flag) force->kspace->setup();
}

/* ---------------------------------------------------------------------- */

void FixNPTAsphere::final_integrate()
{
  int i;
  double dtfm;

  // update v of only atoms in group

  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  if (which == NOBIAS) { 
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	dtfm = dtf / mass[type[i]];
	v[i][0] = (v[i][0] + dtfm*f[i][0]) * factor[0];
	v[i][1] = (v[i][1] + dtfm*f[i][1]) * factor[1];
	v[i][2] = (v[i][2] + dtfm*f[i][2]) * factor[2];
	angmom[i][0] = (angmom[i][0] + dtf * torque[i][0]) * factor_rotate;
	angmom[i][1] = (angmom[i][1] + dtf * torque[i][1]) * factor_rotate;
	angmom[i][2] = (angmom[i][2] + dtf * torque[i][2]) * factor_rotate;
      }
    }
  } else if (which == BIAS) {
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
	temperature->remove_bias(i,v[i]);
	dtfm = dtf / mass[type[i]];
	v[i][0] = (v[i][0] + dtfm*f[i][0]) * factor[0];
	v[i][1] = (v[i][1] + dtfm*f[i][1]) * factor[1];
	v[i][2] = (v[i][2] + dtfm*f[i][2]) * factor[2];
	temperature->restore_bias(i,v[i]);
	angmom[i][0] = (angmom[i][0] + dtf * torque[i][0]) * factor_rotate;
	angmom[i][1] = (angmom[i][1] + dtf * torque[i][1]) * factor_rotate;
	angmom[i][2] = (angmom[i][2] + dtf * torque[i][2]) * factor_rotate;
      }
    }
  }

  // compute new T,P

  t_current = temperature->compute_scalar();
  if (press_couple == 0) {
    double tmp = pressure->compute_scalar();
  } else {
    temperature->compute_vector();
    pressure->compute_vector();
  }
  couple();

  // trigger virial computation on next timestep

  pressure->addstep(update->ntimestep+1);

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;

  // update omega_dot
  // for non-varying dims, p_freq is 0.0, so omega_dot doesn't change

  double f_omega,volume;
  if (dimension == 3) volume = domain->xprd*domain->yprd*domain->zprd;
  else volume = domain->xprd*domain->yprd;
  double denskt = atom->natoms*boltz*t_target / volume * nktv2p;

  for (i = 0; i < 3; i++) {
    f_omega = p_freq[i]*p_freq[i] * (p_current[i]-p_target[i])/denskt;
    omega_dot[i] += f_omega*dthalf;
    omega_dot[i] *= drag_factor;
  }
}

/* ----------------------------------------------------------------------
   Richardson iteration to update quaternion accurately
------------------------------------------------------------------------- */

void FixNPTAsphere::richardson(double *q, double *m, double *moments)
{
  // compute omega at 1/2 step from m at 1/2 step and q at 0

  double w[3];
  omega_from_mq(q,m,moments,w);

  // full update from dq/dt = 1/2 w q

  double wq[4];
  MathExtra::multiply_vec_quat(w,q,wq);

  double qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];
  MathExtra::normalize4(qfull);

  // 1st half of update from dq/dt = 1/2 w q

  double qhalf[4];
  qhalf[0] = q[0] + 0.5*dtq * wq[0];
  qhalf[1] = q[1] + 0.5*dtq * wq[1];
  qhalf[2] = q[2] + 0.5*dtq * wq[2];
  qhalf[3] = q[3] + 0.5*dtq * wq[3];
  MathExtra::normalize4(qhalf);

  // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
  // recompute wq

  omega_from_mq(qhalf,m,moments,w);
  MathExtra::multiply_vec_quat(w,qhalf,wq);

  // 2nd half of update from dq/dt = 1/2 w q

  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];
  MathExtra::normalize4(qhalf);

  // corrected Richardson update

  q[0] = 2.0*qhalf[0] - qfull[0];
  q[1] = 2.0*qhalf[1] - qfull[1];
  q[2] = 2.0*qhalf[2] - qfull[2];
  q[3] = 2.0*qhalf[3] - qfull[3];
  MathExtra::normalize4(q);
}

/* ----------------------------------------------------------------------
   compute omega from angular momentum
   w = omega = angular velocity in space frame
   wbody = angular velocity in body frame
   project space-frame angular momentum onto body axes
     and divide by principal moments
------------------------------------------------------------------------- */

void FixNPTAsphere::omega_from_mq(double *q, double *m, double *inertia,
				  double *w)
{
  double rot[3][3];
  MathExtra::quat_to_mat(q,rot);
  
  double wbody[3];
  MathExtra::transpose_times_column3(rot,m,wbody);
  wbody[0] /= inertia[0];
  wbody[1] /= inertia[1];
  wbody[2] /= inertia[2];
  MathExtra::times_column3(rot,wbody,w);
}


/* ----------------------------------------------------------------------
   principal moments of inertia for ellipsoids
------------------------------------------------------------------------- */

void FixNPTAsphere::calculate_inertia()
{
  double *mass = atom->mass;
  double **shape = atom->shape;
  
  for (int i = 1; i <= atom->ntypes; i++) {
    inertia[i][0] = 0.2*mass[i] *
      (shape[i][1]*shape[i][1]+shape[i][2]*shape[i][2]);
    inertia[i][1] = 0.2*mass[i] *
      (shape[i][0]*shape[i][0]+shape[i][2]*shape[i][2]);
    inertia[i][2] = 0.2*mass[i] * 
      (shape[i][0]*shape[i][0]+shape[i][1]*shape[i][1]);
  }
}
