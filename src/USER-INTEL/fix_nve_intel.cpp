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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "fix_nve_intel.h"
#include "atom.h"
#include "force.h"
#include "intel_preprocess.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEIntel::FixNVEIntel(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  _dtfm = 0;
  _nlocal3 = 0;
  _nlocal_max = 0;
}

/* ---------------------------------------------------------------------- */

FixNVEIntel::~FixNVEIntel()
{
  memory->destroy(_dtfm);
}

/* ---------------------------------------------------------------------- */

void FixNVEIntel::setup(int vflag)
{
  FixNVE::setup(vflag);
  _nlocal3 = 3 * atom->nlocal;
  if (atom->ntypes > 1) reset_dt();
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEIntel::initial_integrate(int /*vflag*/)
{
  // update v and x of atoms in group

  double * _noalias const x = atom->x[0];
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];

  if (igroup == 0 && atom->ntypes == 1 && !atom->rmass) {
    const double dtfm = dtf / atom->mass[1];
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++) {
      v[i] += dtfm * f[i];
      x[i] += dtv * v[i];
    }
  } else if (igroup == 0) {
    if (neighbor->ago == 0) reset_dt();
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++) {
      v[i] += _dtfm[i] * f[i];
      x[i] += dtv * v[i];
    }
  } else {
    if (neighbor->ago == 0) reset_dt();
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++) {
      if (_dtfm[i] != 0.0) {
        v[i] += _dtfm[i] * f[i];
        x[i] += dtv * v[i];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEIntel::final_integrate()
{
  // update v of atoms in group
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];

  if (igroup == 0 && atom->ntypes == 1 && !atom->rmass) {
    _nlocal3 = 3 * atom->nlocal;
    const double dtfm = dtf / atom->mass[1];
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++)
      v[i] += dtfm * f[i];
  } else if (igroup == 0) {
    if (neighbor->ago == 0) reset_dt();
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++) {
      v[i] += _dtfm[i] * f[i];
    }
  } else {
    if (neighbor->ago == 0) reset_dt();
    #if defined(LMP_SIMD_COMPILER)
    #pragma vector aligned
    #pragma simd
    #endif
    for (int i = 0; i < _nlocal3; i++)
      v[i] += _dtfm[i] * f[i];
  }
}

void FixNVEIntel::reset_dt() {
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  const int * const mask = atom->mask;
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;

  if (nlocal > _nlocal_max) {
    if (_nlocal_max) memory->destroy(_dtfm);
    _nlocal_max = static_cast<int>(1.20 * nlocal);
    memory->create(_dtfm, _nlocal_max * 3, "fix_nve_intel:dtfm");
  }

  _nlocal3 = nlocal * 3;

  if (igroup == 0) {
    if (atom->rmass) {
      const double * const rmass = atom->rmass;
      int n = 0;
      for (int i = 0; i < nlocal; i++) {
        _dtfm[n++] = dtf / rmass[i];
        _dtfm[n++] = dtf / rmass[i];
        _dtfm[n++] = dtf / rmass[i];
      }
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      int n = 0;
      for (int i = 0; i < nlocal; i++) {
        _dtfm[n++] = dtf / mass[type[i]];
        _dtfm[n++] = dtf / mass[type[i]];
        _dtfm[n++] = dtf / mass[type[i]];
      }
    }
  } else {
    if (atom->rmass) {
      const double * const rmass = atom->rmass;
      int n = 0;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          _dtfm[n++] = dtf / rmass[i];
          _dtfm[n++] = dtf / rmass[i];
          _dtfm[n++] = dtf / rmass[i];
        } else {
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
        }
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      int n = 0;
      for (int i = 0; i < nlocal; i++)
        if (mask[i] & groupbit) {
          _dtfm[n++] = dtf / mass[type[i]];
          _dtfm[n++] = dtf / mass[type[i]];
          _dtfm[n++] = dtf / mass[type[i]];
        } else {
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
        }
    }
  }
}

double FixNVEIntel::memory_usage()
{
  return FixNVE::memory_usage() + _nlocal_max * 3 * sizeof(double);
}
