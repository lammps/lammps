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
#include "fix_nvt_asphere.h"
#include "math_extra.h"
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

/* ---------------------------------------------------------------------- */

FixNVTASphere::FixNVTASphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVT(lmp, narg, arg)
{
  if (!atom->quat_flag || !atom->angmom_flag || !atom->torque_flag ||
      !atom->avec->shape_type)
    error->all("Fix nvt/asphere requires atom attributes "
	       "quat, angmom, torque, shape");
}

/* ---------------------------------------------------------------------- */

void FixNVTASphere::initial_integrate(int vflag)
{
  double dtfm;

  double delta = update->ntimestep - update->beginstep;
  delta /= update->endstep - update->beginstep;
  t_target = t_start + delta * (t_stop-t_start);

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
  eta += dtv*eta_dot;
  factor = exp(-dthalf*eta_dot);

  // update v and x of only atoms in NVT group

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

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]];
      v[i][0] = v[i][0]*factor + dtfm*f[i][0];
      v[i][1] = v[i][1]*factor + dtfm*f[i][1];
      v[i][2] = v[i][2]*factor + dtfm*f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // update angular momentum by 1/2 step
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion
      
      angmom[i][0] = angmom[i][0] * factor + dtf * torque[i][0];
      angmom[i][1] = angmom[i][1] * factor + dtf * torque[i][1];
      angmom[i][2] = angmom[i][2] * factor + dtf * torque[i][2];

      double inertia[3];
      calculate_inertia(atom->mass[type[i]],atom->shape[type[i]],inertia);
      richardson(quat[i],angmom[i],inertia);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVTASphere::final_integrate()
{
  double dtfm;

  // update v of only atoms in NVT group

  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      dtfm = dtf / mass[type[i]] * factor;
      v[i][0] = v[i][0]*factor + dtfm*f[i][0];
      v[i][1] = v[i][1]*factor + dtfm*f[i][1];
      v[i][2] = v[i][2]*factor + dtfm*f[i][2];

      angmom[i][0] = angmom[i][0] * factor + dtf * torque[i][0];
      angmom[i][1] = angmom[i][1] * factor + dtf * torque[i][1];
      angmom[i][2] = angmom[i][2] * factor + dtf * torque[i][2];
    }
  }

  // compute current T

  t_current = temperature->compute_scalar();

  // update eta_dot

  f_eta = t_freq*t_freq * (t_current/t_target - 1.0);
  eta_dot += f_eta*dthalf;
  eta_dot *= drag_factor;
}

/* ----------------------------------------------------------------------
   Richardson iteration to update quaternion accurately
------------------------------------------------------------------------- */

void FixNVTASphere::richardson(double *q, double *m, double *moments)
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

void FixNVTASphere::omega_from_mq(double *q, double *m, double *inertia,
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
   calculate the moment of inertia for an ELLIPSOID, from mass and radii
   shape = x,y,z radii in body frame
------------------------------------------------------------------------- */

void FixNVTASphere::calculate_inertia(double mass, double *shape,
				      double *inertia)
{
  inertia[0] = mass*(shape[1]*shape[1]+shape[2]*shape[2])/5.0;
  inertia[1] = mass*(shape[0]*shape[0]+shape[2]*shape[2])/5.0;
  inertia[2] = mass*(shape[0]*shape[0]+shape[1]*shape[1])/5.0;
}
