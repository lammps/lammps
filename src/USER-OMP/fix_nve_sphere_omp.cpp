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

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "fix_nve_sphere_omp.h"
#include "atom.h"
#include "atom_vec.h"
#include "update.h"
#include "respa.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.4          // moment of inertia prefactor for sphere

enum{NONE,DIPOLE};

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void FixNVESphereOMP::initial_integrate(int vflag)
{
  double * const * const x = atom->x;
  double * const * const v = atom->v;
  const double * const * const f = atom->f;
  double * const * const omega = atom->omega;
  const double * const * const torque = atom->torque;
  const double * const radius = atom->radius;
  const double * const rmass = atom->rmass;
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  // set timestep here since dt may have changed or come via rRESPA
  const double dtfrotate = dtf / INERTIA;

  // update v,x,omega for all particles
  // d_omega/dt = torque / inertia
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none)
#endif
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];
      x[i][0] += dtv * v[i][0];
      x[i][1] += dtv * v[i][1];
      x[i][2] += dtv * v[i][2];

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
  }

  // update mu for dipoles
  // d_mu/dt = omega cross mu
  // renormalize mu to dipole length

  if (extra == DIPOLE) {
    double * const * const mu = atom->mu;
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none)
#endif
    for (i = 0; i < nlocal; i++) {
      double g0,g1,g2,msq,scale;
      if (mask[i] & groupbit) {
        if (mu[i][3] > 0.0) {
          g0 = mu[i][0] + dtv * (omega[i][1]*mu[i][2]-omega[i][2]*mu[i][1]);
          g1 = mu[i][1] + dtv * (omega[i][2]*mu[i][0]-omega[i][0]*mu[i][2]);
          g2 = mu[i][2] + dtv * (omega[i][0]*mu[i][1]-omega[i][1]*mu[i][0]);
          msq = g0*g0 + g1*g1 + g2*g2;
          scale = mu[i][3]/sqrt(msq);
          mu[i][0] = g0*scale;
          mu[i][1] = g1*scale;
          mu[i][2] = g2*scale;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVESphereOMP::final_integrate()
{
  double * const * const v = atom->v;
  const double * const * const f = atom->f;
  double * const * const omega = atom->omega;
  const double * const * const torque = atom->torque;
  const double * const rmass = atom->rmass;
  const double * const radius = atom->radius;
  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  // set timestep here since dt may have changed or come via rRESPA

  const double dtfrotate = dtf / INERTIA;

  // update v,omega for all particles
  // d_omega/dt = torque / inertia

#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none)
#endif
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i][0] += dtfm * f[i][0];
      v[i][1] += dtfm * f[i][1];
      v[i][2] += dtfm * f[i][2];

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i][0] += dtirotate * torque[i][0];
      omega[i][1] += dtirotate * torque[i][1];
      omega[i][2] += dtirotate * torque[i][2];
    }
}
