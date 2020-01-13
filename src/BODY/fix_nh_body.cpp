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
   Contributing author: Trung Dac Nguyen (ndactrung@gmail.com)
   based on FixNHAsphere
------------------------------------------------------------------------- */

#include "fix_nh_body.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_body.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNHBody::FixNHBody(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void FixNHBody::init()
{
  avec = (AtomVecBody *) atom->style_match("body");
  if (!avec)
    error->all(FLERR,
               "Compute nvt/nph/npt body requires atom style body");

  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  int *body = atom->body;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (body[i] < 0)
        error->one(FLERR,"Fix nvt/nph/npt body requires bodies");

  FixNH::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of angular momentum
-----------------------------------------------------------------------*/

void FixNHBody::nve_v()
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

void FixNHBody::nve_x()
{
  double omega[3];
  double *quat,*inertia;

  // standard nve_x position update

  FixNH::nve_x();

  AtomVecBody::Bonus *bonus = avec->bonus;
  int *body = atom->body;
  double **angmom = atom->angmom;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  // update quaternion a full step via Richardson iteration

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion

      inertia = bonus[body[i]].inertia;
      quat = bonus[body[i]].quat;
      MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
      MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);
    }
}

/* ----------------------------------------------------------------------
   perform half-step temperature scaling of angular momentum
-----------------------------------------------------------------------*/

void FixNHBody::nh_v_temp()
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
