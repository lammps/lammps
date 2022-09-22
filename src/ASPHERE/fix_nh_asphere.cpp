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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#include "fix_nh_asphere.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNHAsphere::FixNHAsphere(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void FixNHAsphere::init()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  if (!avec)
    error->all(FLERR,
               "Compute nvt/nph/npt asphere requires atom style ellipsoid");

  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nvt/nph/npt asphere requires extended particles");

  FixNH::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of angular momentum
-----------------------------------------------------------------------*/

void FixNHAsphere::nve_v()
{
  // standard nve_v velocity update

  FixNH::nve_v();

  double **angmom = atom->angmom;
  double **torque = atom->torque;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // update angular momentum by 1/2 step for all particles

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      angmom[i][0] += dtf*torque[i][0];
      angmom[i][1] += dtf*torque[i][1];
      angmom[i][2] += dtf*torque[i][2];
    }
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of orientation
-----------------------------------------------------------------------*/

void FixNHAsphere::nve_x()
{
  double omega[3];

  // standard nve_x position update

  FixNH::nve_x();

  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double **angmom = atom->angmom;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  // update quaternion a full step via Richardson iteration
  // returns new normalized quaternion
  // principal moments of inertia

  double *shape,*quat;
  double inertia[3];

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      // principal moments of inertia

      shape = bonus[ellipsoid[i]].shape;
      quat = bonus[ellipsoid[i]].quat;

      inertia[0] = rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]) / 5.0;
      inertia[1] = rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]) / 5.0;
      inertia[2] = rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]) / 5.0;

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion

      MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
      MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);
    }
}

/* ----------------------------------------------------------------------
   perform half-step temperature scaling of angular momentum
-----------------------------------------------------------------------*/

void FixNHAsphere::nh_v_temp()
{
  // standard nh_v_temp scaling

  FixNH::nh_v_temp();

  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      angmom[i][0] *= factor_eta;
      angmom[i][1] *= factor_eta;
      angmom[i][2] *= factor_eta;
    }
  }
}
