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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "fix_nve_asphere_intel.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "error.h"
#include "force.h"
#include "math_extra_intel.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

#define INERTIA 0.2          // moment of inertia prefactor for ellipsoid

/* ---------------------------------------------------------------------- */

FixNVEAsphereIntel::FixNVEAsphereIntel(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  _dtfm = nullptr;
  _nlocal3 = 0;
  _nlocal_max = 0;
  _inertia0 = nullptr;
  _inertia1 = nullptr;
  _inertia2 = nullptr;
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereIntel::init()
{
  avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
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

void FixNVEAsphereIntel::setup(int vflag)
{
  FixNVE::setup(vflag);
  reset_dt();
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereIntel::initial_integrate(int /*vflag*/)
{
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;
  double * _noalias const x = atom->x[0];
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];
  int *mask = atom->mask;

  double **angmom = atom->angmom;
  double **torque = atom->torque;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  // set timestep here since dt may have changed or come via rRESPA

  dtq = 0.5 * dtv;

  #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
  #pragma omp simd
#else
  #pragma simd
#endif
  #pragma vector aligned
  #endif
  for (int i = 0; i < _nlocal3; i++) {
    v[i] += _dtfm[i] * f[i];
    x[i] += dtv * v[i];
  }

  // update angular momentum by 1/2 step
  if (igroup == 0) {
    #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
    #pragma omp simd
#else
    #pragma simd
#endif
    #pragma vector aligned
    #endif
    for (int i = 0; i < nlocal; i++) {
      double *quat = bonus[ellipsoid[i]].quat;
      ME_omega_richardson(dtf, dtq, angmom[i], quat, torque[i], _inertia0[i],
                          _inertia1[i], _inertia2[i]);
    }
  } else {
    #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
    #pragma omp simd
#else
    #pragma simd
#endif
    #pragma vector aligned
    #endif
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        double *quat = bonus[ellipsoid[i]].quat;
        ME_omega_richardson(dtf, dtq, angmom[i], quat, torque[i], _inertia0[i],
                            _inertia1[i], _inertia2[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEAsphereIntel::final_integrate()
{
  if (neighbor->ago == 0) reset_dt();

  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];
  double * _noalias const angmom = atom->angmom[0];
  const double * _noalias const torque = atom->torque[0];

  #if defined(LMP_SIMD_COMPILER)
#if defined(USE_OMP_SIMD)
  #pragma omp simd
#else
  #pragma simd
#endif
  #pragma vector aligned
  #endif
  for (int i = 0; i < _nlocal3; i++) {
    v[i] += _dtfm[i] * f[i];
    angmom[i] += dtf * torque[i];
  }
}

void FixNVEAsphereIntel::reset_dt() {
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;

  if (nlocal > _nlocal_max) {
    if (_nlocal_max) {
      memory->destroy(_dtfm);
      memory->destroy(_inertia0);
      memory->destroy(_inertia1);
      memory->destroy(_inertia2);
    }
    _nlocal_max = static_cast<int>(1.20 * nlocal);
    memory->create(_dtfm, _nlocal_max * 3, "fix_nve_intel:dtfm");
    memory->create(_inertia0, _nlocal_max * 3, "fix_nve_intel:inertia0");
    memory->create(_inertia1, _nlocal_max * 3, "fix_nve_intel:inertia1");
    memory->create(_inertia2, _nlocal_max * 3, "fix_nve_intel:inertia2");
  }

  _nlocal3 = nlocal * 3;

  if (igroup == 0) {
    const double * const rmass = atom->rmass;
    int n = 0;
    for (int i = 0; i < nlocal; i++) {
      _dtfm[n++] = dtf / rmass[i];
      _dtfm[n++] = dtf / rmass[i];
      _dtfm[n++] = dtf / rmass[i];
      double *shape = bonus[ellipsoid[i]].shape;
      double idot = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
      if (idot != 0.0) idot = 1.0 / idot;
      _inertia0[i] = idot;
      idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
      if (idot != 0.0) idot = 1.0 / idot;
      _inertia1[i] = idot;
      idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
      if (idot != 0.0) idot = 1.0 / idot;
      _inertia2[i] = idot;
    }
  } else {
    const double * const rmass = atom->rmass;
    int n = 0;
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        _dtfm[n++] = dtf / rmass[i];
        _dtfm[n++] = dtf / rmass[i];
        _dtfm[n++] = dtf / rmass[i];
        double *shape = bonus[ellipsoid[i]].shape;
        double idot = INERTIA*rmass[i] * (shape[1]*shape[1]+shape[2]*shape[2]);
        if (idot != 0.0) idot = 1.0 / idot;
        _inertia0[i] = idot;
        idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[2]*shape[2]);
        if (idot != 0.0) idot = 1.0 / idot;
        _inertia1[i] = idot;
        idot = INERTIA*rmass[i] * (shape[0]*shape[0]+shape[1]*shape[1]);
        if (idot != 0.0) idot = 1.0 / idot;
        _inertia2[i] = idot;
      } else {
        _dtfm[n++] = 0.0;
        _dtfm[n++] = 0.0;
        _dtfm[n++] = 0.0;
      }
    }
  }
}
double FixNVEAsphereIntel::memory_usage()
{
  return FixNVE::memory_usage() + _nlocal_max * 12 * sizeof(double);
}

