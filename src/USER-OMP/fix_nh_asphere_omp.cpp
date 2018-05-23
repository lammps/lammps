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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <cstring>
#include <cstdlib>
#include <cmath>
#include "math_extra.h"
#include "fix_nh_asphere_omp.h"
#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "compute.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};

typedef struct { double x,y,z; } dbl3_t;

/* ---------------------------------------------------------------------- */

FixNHAsphereOMP::FixNHAsphereOMP(LAMMPS *lmp, int narg, char **arg) :
  FixNHOMP(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

void FixNHAsphereOMP::init()
{
  avec = (AtomVecEllipsoid *) atom->style_match("ellipsoid");
  if (!avec)
    error->all(FLERR,"Compute nvt/nph/npt asphere requires atom style ellipsoid");

  // check that all particles are finite-size
  // no point particles allowed, spherical is OK

  int *ellipsoid = atom->ellipsoid;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (ellipsoid[i] < 0)
        error->one(FLERR,"Fix nvt/nph/npt asphere requires extended particles");

  FixNHOMP::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of angular momentum and COM velocity
-----------------------------------------------------------------------*/

void FixNHAsphereOMP::nve_v()
{
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  dbl3_t * _noalias const angmom = (dbl3_t *) atom->angmom[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const dbl3_t * _noalias const torque = (dbl3_t *) atom->torque[0];
  const double * _noalias const rmass = atom->rmass;
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  // standard nve_v velocity update. for efficiency the loop is
  // merged with FixNHOMP instead of calling it for the COM update.

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i].x += dtfm*f[i].x;
      v[i].y += dtfm*f[i].y;
      v[i].z += dtfm*f[i].z;
      angmom[i].x += dtf*torque[i].x;
      angmom[i].y += dtf*torque[i].y;
      angmom[i].z += dtf*torque[i].z;
    }
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of position and orientation
-----------------------------------------------------------------------*/

void FixNHAsphereOMP::nve_x()
{
  dbl3_t * _noalias const x = (dbl3_t *) atom->x[0];
  const dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  dbl3_t * _noalias const angmom = (dbl3_t *) atom->angmom[0];
  const double * _noalias const rmass = atom->rmass;
  const int * _noalias const mask = atom->mask;
  AtomVecEllipsoid::Bonus * _noalias const bonus = avec->bonus;
  const int * _noalias const ellipsoid = atom->ellipsoid;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  // standard nve_x position update and
  // update quaternion a full step via Richardson iteration
  // returns new normalized quaternion
  // principal moments of inertia

#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      double omega[3], inertia[3];

      x[i].x += dtv * v[i].x;
      x[i].y += dtv * v[i].y;
      x[i].z += dtv * v[i].z;

      // principal moments of inertia

      const double * const shape = bonus[ellipsoid[i]].shape;
      double * const quat = bonus[ellipsoid[i]].quat;

      inertia[0] = rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]) / 5.0;
      inertia[1] = rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]) / 5.0;
      inertia[2] = rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]) / 5.0;

      // compute omega at 1/2 step from angmom at 1/2 step and current q
      // update quaternion a full step via Richardson iteration
      // returns new normalized quaternion

      MathExtra::mq_to_omega(&angmom[i].x,quat,inertia,omega);
      MathExtra::richardson(quat,&angmom[i].x,omega,inertia,dtq);
    }
}

/* ----------------------------------------------------------------------
   perform half-step temperature scaling of angular momentum
-----------------------------------------------------------------------*/

void FixNHAsphereOMP::nh_v_temp()
{
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  dbl3_t * _noalias const angmom = (dbl3_t *) atom->angmom[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;
  int i;

  if (which == NOBIAS) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i].x *= factor_eta;
        v[i].y *= factor_eta;
        v[i].z *= factor_eta;
        angmom[i].x *= factor_eta;
        angmom[i].y *= factor_eta;
        angmom[i].z *= factor_eta;
      }
    }
  } else if (which == BIAS) {
#if defined(_OPENMP)
#pragma omp parallel for default(none) private(i) schedule(static)
#endif
    for (i = 0; i < nlocal; i++) {
      double buf[3];
      if (mask[i] & groupbit) {
        temperature->remove_bias_thr(i,&v[i].x,buf);
        v[i].x *= factor_eta;
        v[i].y *= factor_eta;
        v[i].z *= factor_eta;
        temperature->restore_bias_thr(i,&v[i].x,buf);
        angmom[i].x *= factor_eta;
        angmom[i].y *= factor_eta;
        angmom[i].z *= factor_eta;
      }
    }
  }
}
