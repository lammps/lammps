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
#include "stdio.h"
#include "string.h"
#include "fix_nve_asphere.h"
#include "math_extra.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "force.h"
#include "update.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixNVEAsphere::FixNVEAsphere(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg) {}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::init()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Compute nve/asphere requires atom style ellipsoid");

  // check that all particles are finite-size ellipsoids
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nve/asphere requires extended particles");

  FixNVE::init();
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::initial_integrate(int vflag)
{
  double dtfm;
  double inertia[3],omega[3];
  double *shape,*quat;

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

  dtq = 0.5 * dtv;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      // update angular momentum by 1/2 step

      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];

      // principal moments of inertia

      shape = bonus[ellipsoid[i]].shape;
      quat = bonus[ellipsoid[i]].quat;

      inertia[0] = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      inertia[1] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      inertia[2] = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion

      MathExtra::mq_to_omega(angmom[i],quat,inertia,omega);
      MathExtra::richardson(quat,angmom[i],omega,inertia,dtq);
    }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphere::final_integrate()
{
  double dtfm;

  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      angmom[i][0] += dtf * torque[i][0];
      angmom[i][1] += dtf * torque[i][1];
      angmom[i][2] += dtf * torque[i][2];
    }
}
