// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "fix_nve_dotc_langevin.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "math_extra.h"
#include "random_mars.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixNVEDotcLangevin::FixNVEDotcLangevin(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  if (narg != 9) error->all(FLERR,"Illegal fix nve/dotc/langevin command");

  t_start = utils::numeric(FLERR,arg[3],false,lmp);
  t_target = t_start;
  t_stop = utils::numeric(FLERR,arg[4],false,lmp);
  t_period = utils::numeric(FLERR,arg[5],false,lmp);
  if (t_period <= 0.0) error->all(FLERR,"Fix nve/dotc/langevin period must be > 0.0");
  gamma = 1.0/t_period;
  seed = utils::inumeric(FLERR,arg[6],false,lmp);
  if (seed <= 0) error->all(FLERR,"Illegal fix nve/dotc/langevin command");

  if (strcmp(arg[7],"angmom") == 0) {
    if (9 > narg) error->all(FLERR,"Illegal fix nve/dotc/langevin command");
    if (strcmp(arg[8],"no") == 0) {
      ascale = 0.0;
      Gamma = 0.0;
    }
    else {
      ascale = utils::numeric(FLERR,arg[8],false,lmp);
      Gamma = gamma * ascale;
    }

  }

  // initialize Marsaglia RNG with processor-unique seed

  random = new RanMars(lmp,seed + comm->me);

}

/* ---------------------------------------------------------------------- */

FixNVEDotcLangevin::~FixNVEDotcLangevin()
{

  delete random;

}


/* ---------------------------------------------------------------------- */

void FixNVEDotcLangevin::init()
{

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");

  if (!avec)
    error->all(FLERR,"Fix nve/dotc/langevin requires atom style ellipsoid");

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nve/dotc/langevin requires extended particles");


  // set prefactor
  gfactor1 = exp(-gamma*update->dt);

  // set square root of temperature
  compute_target();

  FixNVE::init();
}

/* ----------------------------------------------------------------------
   set current t_target and t_sqrt
------------------------------------------------------------------------- */

void FixNVEDotcLangevin::compute_target()
{
  double delta = update->ntimestep - update->beginstep;
  if (delta != 0.0) delta /= update->endstep - update->beginstep;

  // Only homogeneous temperature supported
  t_target = t_start + delta * (t_stop-t_start);
  tsqrt = sqrt(t_target);

}


/* ---------------------------------------------------------------------- */

void FixNVEDotcLangevin::initial_integrate(int /*vflag*/)
{
  double *shape,*quat;
  double fquat[4],conjqm[4],inertia[3];
  double slq_conjqm[3];

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dt = update->dt;
  dthlf = 0.5 * dt;
  dtqrt = 0.25 * dt;

  // set square root of temperature
  compute_target();

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      dthlfm = dthlf / rmass[i];
      quat = bonus[ellipsoid[i]].quat;
      shape = bonus[ellipsoid[i]].shape;

      // update momentum by 1/2 step
      v[i][0] += dthlfm * f[i][0];
      v[i][1] += dthlfm * f[i][1];
      v[i][2] += dthlfm * f[i][2];

      // update position by 1/2 step
      x[i][0] += dthlf * v[i][0];
      x[i][1] += dthlf * v[i][1];
      x[i][2] += dthlf * v[i][2];

      // convert angular momentum and torque in space frame into
      // quaternion 4-momentum and 1/2 of 4-torque in body frame
      vec3_to_vec4(quat,angmom[i],conjqm);
      conjqm[0] *= 2.0;
      conjqm[1] *= 2.0;
      conjqm[2] *= 2.0;
      conjqm[3] *= 2.0;
      vec3_to_vec4(quat,torque[i],fquat);

      // update quaternion 4-momentum by 1/2 step
      conjqm[0] += dt * fquat[0];
      conjqm[1] += dt * fquat[1];
      conjqm[2] += dt * fquat[2];
      conjqm[3] += dt * fquat[3];

      // principal moments of inertia
      inertia[0] = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      inertia[1] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      inertia[2] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);

      M = inertia[0]*inertia[1]*inertia[2];
      M /= inertia[1]*inertia[2]+inertia[0]*inertia[2]+inertia[0]*inertia[1];

      // set prefactors
      // factors 12 and 48 reflect the variance of the uniform distribution:
      // var = 1/12*(b-a)^2
      gfactor2 = sqrt(12.0*(1.0-gfactor1*gfactor1)/rmass[i])*tsqrt;

      gfactor3[0] = exp(-Gamma*M*dt/inertia[0]);
      gfactor3[1] = exp(-Gamma*M*dt/inertia[1]);
      gfactor3[2] = exp(-Gamma*M*dt/inertia[2]);

      gfactor4[0] = sqrt(48.0*inertia[0]*(1.0-gfactor3[0]*gfactor3[0]))*tsqrt;
      gfactor4[1] = sqrt(48.0*inertia[1]*(1.0-gfactor3[1]*gfactor3[1]))*tsqrt;
      gfactor4[2] = sqrt(48.0*inertia[2]*(1.0-gfactor3[2]*gfactor3[2]))*tsqrt;

      // rotate quaternion and quaternion 4-momentum by 1/2 step
      no_squish_rotate(3,conjqm,quat,inertia,dtqrt);
      no_squish_rotate(2,conjqm,quat,inertia,dtqrt);
      no_squish_rotate(1,conjqm,quat,inertia,dthlf);
      no_squish_rotate(2,conjqm,quat,inertia,dtqrt);
      no_squish_rotate(3,conjqm,quat,inertia,dtqrt);

      // apply stochastic force to velocities
      v[i][0] = v[i][0] * gfactor1 + gfactor2 * (random->uniform()-0.5);
      v[i][1] = v[i][1] * gfactor1 + gfactor2 * (random->uniform()-0.5);
      v[i][2] = v[i][2] * gfactor1 + gfactor2 * (random->uniform()-0.5);

      // update position by 1/2 step
      x[i][0] += dthlf * v[i][0];
      x[i][1] += dthlf * v[i][1];
      x[i][2] += dthlf * v[i][2];

      // apply stochastic force to quaternion 4-momentum
      slq_conjqm[0] = -quat[1]*conjqm[0] + quat[0]*conjqm[1] + quat[3]*conjqm[2] - quat[2]*conjqm[3];
      slq_conjqm[1] = -quat[2]*conjqm[0] - quat[3]*conjqm[1] + quat[0]*conjqm[2] + quat[1]*conjqm[3];
      slq_conjqm[2] = -quat[3]*conjqm[0] + quat[2]*conjqm[1] - quat[1]*conjqm[2] + quat[0]*conjqm[3];

      gfactor5[0] = gfactor3[0] * slq_conjqm[0] + gfactor4[0] * (random->uniform()-0.5);
      gfactor5[1] = gfactor3[1] * slq_conjqm[1] + gfactor4[1] * (random->uniform()-0.5);
      gfactor5[2] = gfactor3[2] * slq_conjqm[2] + gfactor4[2] * (random->uniform()-0.5);

      conjqm[0] = -quat[1] * gfactor5[0] - quat[2] * gfactor5[1] - quat[3] * gfactor5[2];
      conjqm[1] =  quat[0] * gfactor5[0] - quat[3] * gfactor5[1] + quat[2] * gfactor5[2];
      conjqm[2] =  quat[3] * gfactor5[0] + quat[0] * gfactor5[1] - quat[1] * gfactor5[2];
      conjqm[3] = -quat[2] * gfactor5[0] + quat[1] * gfactor5[1] + quat[0] * gfactor5[2];

      // rotate quaternion and quaternion 4-momentum by 1/2 step
      no_squish_rotate(3,conjqm,quat,inertia,dtqrt);
      no_squish_rotate(2,conjqm,quat,inertia,dtqrt);
      no_squish_rotate(1,conjqm,quat,inertia,dthlf);
      no_squish_rotate(2,conjqm,quat,inertia,dtqrt);
      no_squish_rotate(3,conjqm,quat,inertia,dtqrt);
      qnormalize(quat);

      // convert quaternion 4-momentum in body frame back to angular momentum in space frame
      vec4_to_vec3(quat,conjqm,angmom[i]);

      angmom[i][0] *= 0.5;
      angmom[i][1] *= 0.5;
      angmom[i][2] *= 0.5;

    }

}

/* ---------------------------------------------------------------------- */

void FixNVEDotcLangevin::final_integrate()
{

  double *quat;
  double fquat[4],conjqm[4];
  double conjqm_dot_quat;

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dt = update->dt;
  dthlf = 0.5 * dt;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      dthlfm = dthlf / rmass[i];
      quat = bonus[ellipsoid[i]].quat;

      // update momentum by 1/2 step
      v[i][0] += dthlfm * f[i][0];
      v[i][1] += dthlfm * f[i][1];
      v[i][2] += dthlfm * f[i][2];

      // convert angular momentum and torque in space frame into
      // quaternion 4-momentum and 1/2 of 4-torque in body frame
      vec3_to_vec4(quat,angmom[i],conjqm);
      conjqm[0] *= 2.0;
      conjqm[1] *= 2.0;
      conjqm[2] *= 2.0;
      conjqm[3] *= 2.0;
      vec3_to_vec4(quat,torque[i],fquat);

      // update quaternion 4-momentum by 1/2 step
      conjqm[0] += dt * fquat[0];
      conjqm[1] += dt * fquat[1];
      conjqm[2] += dt * fquat[2];
      conjqm[3] += dt * fquat[3];

      // subtract component parallel to quaternion for improved numerical accuracy

      conjqm_dot_quat = conjqm[0]*quat[0] + conjqm[1]*quat[1] + conjqm[2]*quat[2] + conjqm[3]*quat[3];

      conjqm[0] -= conjqm_dot_quat * quat[0];
      conjqm[1] -= conjqm_dot_quat * quat[1];
      conjqm[2] -= conjqm_dot_quat * quat[2];
      conjqm[3] -= conjqm_dot_quat * quat[3];

      // convert quaternion 4-momentum in body frame back to angular momentum in space frame
      vec4_to_vec3(quat,conjqm,angmom[i]);

      angmom[i][0] *= 0.5;
      angmom[i][1] *= 0.5;
      angmom[i][2] *= 0.5;

    }
}
