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
#include "math_extra.h"
#include "fix_nh_asphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

enum{NOBIAS,BIAS};

/* ---------------------------------------------------------------------- */

FixNHAsphere::FixNHAsphere(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  inertia = 
    memory->create_2d_double_array(atom->ntypes+1,3,"fix_nvt_asphere:inertia");

  if (!atom->quat_flag || !atom->angmom_flag || !atom->torque_flag ||
      !atom->avec->shape_type)
    error->all("Fix nvt/nph/npt asphere requires atom attributes "
	       "quat, angmom, torque, shape");
  if (atom->radius_flag || atom->rmass_flag)
    error->all("Fix nvt/nph/npt asphere cannot be used with atom attributes "
	       "diameter or rmass");
}

/* ---------------------------------------------------------------------- */

FixNHAsphere::~FixNHAsphere()
{
  memory->destroy_2d_double_array(inertia);
}

/* ---------------------------------------------------------------------- */

void FixNHAsphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  double **shape = atom->shape;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (shape[type[i]][0] == 0.0)
	error->one("Fix nvt/nph/npt asphere requires extended particles");

  FixNH::init();
  calculate_inertia();
}

/* ----------------------------------------------------------------------
   Richardson iteration to update quaternion accurately
------------------------------------------------------------------------- */

void FixNHAsphere::richardson(double *q, double *m, double *moments)
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

void FixNHAsphere::omega_from_mq(double *q, double *m, double *inertia,
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

void FixNHAsphere::calculate_inertia()
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

/* ----------------------------------------------------------------------
   perform half-step update of velocities 
-----------------------------------------------------------------------*/

void FixNHAsphere::nve_v()
{
  // standard nhc_nve_v velocity update

  FixNH::nve_v();

  int *type = atom->type;
  double **quat = atom->quat;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  // update angular momentum by 1/2 step for all particles
  // update quaternion a full step via Richardson iteration
  // returns new normalized quaternion

  for (int i = 0; i < nlocal; i++) {    
    if (mask[i] & groupbit) {
      angmom[i][0] += dtf*torque[i][0];
      angmom[i][1] += dtf*torque[i][1];
      angmom[i][2] += dtf*torque[i][2];
		
      richardson(quat[i],angmom[i],inertia[type[i]]);
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step temperature scaling of velocities
-----------------------------------------------------------------------*/

void FixNHAsphere::nh_v_temp()
{
  // standard nhc_nh_v velocity update

  FixNH::nh_v_temp();

  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double factor_rotate = exp(-dthalf*eta_dot[0]);

  for (int i = 0; i < nlocal; i++) {    
    if (mask[i] & groupbit) {
      angmom[i][0] *= factor_rotate;
      angmom[i][1] *= factor_rotate;
      angmom[i][2] *= factor_rotate;
    }
  }
}
