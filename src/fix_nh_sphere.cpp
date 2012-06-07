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
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

/* ---------------------------------------------------------------------- */

FixNHSphere::FixNHSphere(LAMMPS *lmp, int narg, char **arg) :
  FixNH(lmp, narg, arg)
{
  if (!atom->sphere_flag)
    error->all(FLERR,"Fix nvt/nph/npt sphere requires atom style sphere");
}

/* ---------------------------------------------------------------------- */

void FixNHSphere::init()
{
  // check that all particles are finite-size
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nvt/sphere requires extended particles");

  FixNH::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of rotational velocities
-----------------------------------------------------------------------*/

void FixNHSphere::nve_v()
{
  // standard nve_v velocity update

  FixNH::nve_v();

  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  double dtfrotate = dtf / INERTIA;
  double dtirotate;

  // update omega for all particles
  // d_omega/dt = torque / inertia
  // 4 cases depending on radius vs shape and rmass vs mass

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate*torque[i][0];
      omega[i][1] += dtirotate*torque[i][1];
      omega[i][2] += dtirotate*torque[i][2];
    }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of rotatonal velocities
-----------------------------------------------------------------------*/

void FixNHSphere::nh_v_temp()
{
  // standard nh_v_temp scaling

  FixNH::nh_v_temp();

  double **omega = atom->omega;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      omega[i][0] *= factor_eta;
      omega[i][1] *= factor_eta;
      omega[i][2] *= factor_eta;
    }
  }
}
