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

#include "omp_compat.h"
#include "fix_nh_sphere_omp.h"
#include "atom.h"
#include "compute.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NOBIAS,BIAS};

#define INERTIA 0.4          // moment of inertia prefactor for sphere

typedef struct { double x,y,z; } dbl3_t;

/* ---------------------------------------------------------------------- */

FixNHSphereOMP::FixNHSphereOMP(LAMMPS *lmp, int narg, char **arg) :
  FixNHOMP(lmp, narg, arg)
{
  if (!atom->sphere_flag)
    error->all(FLERR,"Fix nvt/nph/npt sphere requires atom style sphere");
}

/* ---------------------------------------------------------------------- */

void FixNHSphereOMP::init()
{
  // check that all particles are finite-size
  // no point particles allowed

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0)
        error->one(FLERR,"Fix nvt/npt/nph/sphere/omp require extended particles");

  FixNHOMP::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of rotational and COM velocities
-----------------------------------------------------------------------*/

void FixNHSphereOMP::nve_v()
{
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  dbl3_t * _noalias const omega = (dbl3_t *) atom->omega[0];
  const dbl3_t * _noalias const f = (dbl3_t *) atom->f[0];
  const dbl3_t * _noalias const torque = (dbl3_t *) atom->torque[0];
  const double * _noalias const radius = atom->radius;
  const double * _noalias const rmass = atom->rmass;
  const int * _noalias const mask = atom->mask;

  // set timestep here since dt may have changed or come via rRESPA

  const double dtfrotate = dtf / INERTIA;

  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  // standard nve_v velocity update. for efficiency the loop is
  // merged with FixNHOMP instead of calling it for the COM update.

  // update omega for all particles
  // d_omega/dt = torque / inertia
  // 4 cases depending on radius vs shape and rmass vs mass

#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      const double dtfm = dtf / rmass[i];
      v[i].x += dtfm*f[i].x;
      v[i].y += dtfm*f[i].y;
      v[i].z += dtfm*f[i].z;

      const double dtirotate = dtfrotate / (radius[i]*radius[i]*rmass[i]);
      omega[i].x += dtirotate*torque[i].x;
      omega[i].y += dtirotate*torque[i].y;
      omega[i].z += dtirotate*torque[i].z;
    }
  }
}

/* ----------------------------------------------------------------------
   perform half-step scaling of rotatonal velocities
-----------------------------------------------------------------------*/

void FixNHSphereOMP::nh_v_temp()
{
  dbl3_t * _noalias const v = (dbl3_t *) atom->v[0];
  dbl3_t * _noalias const omega = (dbl3_t *) atom->omega[0];
  const int * _noalias const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst : atom->nlocal;

  if (which == NOBIAS) {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        v[i].x *= factor_eta;
        v[i].y *= factor_eta;
        v[i].z *= factor_eta;
        omega[i].x *= factor_eta;
        omega[i].y *= factor_eta;
        omega[i].z *= factor_eta;
      }
    }
  } else if (which == BIAS) {
#if defined(_OPENMP)
#pragma omp parallel for LMP_DEFAULT_NONE schedule(static)
#endif
    for (int i = 0; i < nlocal; i++) {
      double buf[3];
      if (mask[i] & groupbit) {
        temperature->remove_bias_thr(i,&v[i].x,buf);
        v[i].x *= factor_eta;
        v[i].y *= factor_eta;
        v[i].z *= factor_eta;
        temperature->restore_bias_thr(i,&v[i].x,buf);
        omega[i].x *= factor_eta;
        omega[i].y *= factor_eta;
        omega[i].z *= factor_eta;
      }
    }
  }
}
