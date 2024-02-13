// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "fix_nve_dot.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "update.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

static constexpr double INERTIA = 0.2;          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixNVEDot::FixNVEDot(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixNVEDot::init()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec)
    error->all(FLERR,"Compute nve/dot requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nve/dot requires extended particles");

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

void FixNVEDot::initial_integrate(int /*vflag*/)
{
  double *shape,*quat;
  double fquat[4],conjqm[4],inertia[3];

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

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      dthlfm = dthlf / rmass[i];
      quat = bonus[ellipsoid[i]].quat;
      shape = bonus[ellipsoid[i]].shape;

      // update momentum by 1/2 step
      v[i][0] += dthlfm * f[i][0];
      v[i][1] += dthlfm * f[i][1];
      v[i][2] += dthlfm * f[i][2];

      // update position by full step
      x[i][0] += dt * v[i][0];
      x[i][1] += dt * v[i][1];
      x[i][2] += dt * v[i][2];

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

      // rotate quaternion and quaternion 4-momentum by full step
      no_squish_rotate(3,conjqm,quat,inertia,dthlf);
      no_squish_rotate(2,conjqm,quat,inertia,dthlf);
      no_squish_rotate(1,conjqm,quat,inertia,dt);
      no_squish_rotate(2,conjqm,quat,inertia,dthlf);
      no_squish_rotate(3,conjqm,quat,inertia,dthlf);

      qnormalize(quat);

      // convert quaternion 4-momentum in body frame back to angular momentum in space frame
      vec4_to_vec3(quat,conjqm,angmom[i]);

      angmom[i][0] *= 0.5;
      angmom[i][1] *= 0.5;
      angmom[i][2] *= 0.5;

    }
}

/* ---------------------------------------------------------------------- */

void FixNVEDot::final_integrate()
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

      // update momentum
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
