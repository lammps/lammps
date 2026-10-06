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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "fix_nve_gpu.h"

#include "atom.h"
#include "force.h"
#include "gpu_extra.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#if (LAL_USE_OMP == 1)
#include <omp.h>
#endif

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVEGPU::FixNVEGPU(LAMMPS *lmp, int narg, char **arg) :
  FixNVE(lmp, narg, arg)
{
  _dtfm = nullptr;
  _nlocal_max = 0;
}

/* ---------------------------------------------------------------------- */

FixNVEGPU::~FixNVEGPU()
{
  memory->destroy(_dtfm);
}

/* ---------------------------------------------------------------------- */

void FixNVEGPU::setup(int vflag)
{
  FixNVE::setup(vflag);
  if (utils::strmatch(update->integrate_style,"^respa"))
    _respa_on = 1;
  else
    _respa_on = 0;

  // ensure that _dtfm array is initialized if the group is not "all"
  // or there is more than one atom type as that re-ordeted array is used for
  // per-type/per-atom masses and group membership detection.
  if ((igroup != 0) || (atom->ntypes > 1)) reset_dt();
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVEGPU::initial_integrate(int vflag)
{
  if (_respa_on) { FixNVE::initial_integrate(vflag); return; }

  // update v and x of atoms in group

  double * _noalias const x = atom->x[0];
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;

  #if (LAL_USE_OMP == 1)
  #pragma omp parallel
  #endif
  {
    #if (LAL_USE_OMP == 1)
    const int nthreads = omp_get_num_threads();
    const int idelta = nlocal / nthreads + 1;
    const int ifrom = omp_get_thread_num() * idelta;
    const int ito = MIN(ifrom + idelta, nlocal);
    #else
    const int ifrom = 0;
    const int ito = nlocal;
    #endif
    const int ifrom3 = ifrom * 3;
    const int ito3 = ito * 3;

    if (igroup == 0 && atom->ntypes == 1 && !atom->rmass) {
      const double dtfm = dtf / atom->mass[1];
      #if (LAL_USE_OMP_SIMD == 1)
      #pragma omp simd
      #endif
      for (int i = ifrom3; i < ito3; i++) {
        v[i] += dtfm * f[i];
        x[i] += dtv * v[i];
      }
    } else if (igroup == 0) {
      for (int i = ifrom; i < ito; i++) {
        const double dtfm = _dtfm[i];
        const int n = i * 3;
        v[n] += dtfm * f[n];
        v[n+1] += dtfm * f[n+1];
        v[n+2] += dtfm * f[n+2];
        x[n] += dtv * v[n];
        x[n+1] += dtv * v[n+1];
        x[n+2] += dtv * v[n+2];
      }
    } else {
      for (int i = ifrom; i < ito; i++) {
        const double dtfm = _dtfm[i];
        if (dtfm == 0.0) continue;
        const int n = i * 3;
        v[n] += dtfm * f[n];
        v[n+1] += dtfm * f[n+1];
        v[n+2] += dtfm * f[n+2];
        x[n] += dtv * v[n];
        x[n+1] += dtv * v[n+1];
        x[n+2] += dtv * v[n+2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixNVEGPU::final_integrate()
{
  if (_respa_on) { FixNVE::final_integrate(); return; }
  // update v of atoms in group
  double * _noalias const v = atom->v[0];
  const double * _noalias const f = atom->f[0];
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;

  if (neighbor->ago == 0) {
    if (igroup != 0 || atom->ntypes != 1 || atom->rmass) {
      if (nlocal > _nlocal_max) {
        if (_nlocal_max) memory->destroy(_dtfm);
        _nlocal_max = static_cast<int>(1.20 * nlocal);
        memory->create(_dtfm, _nlocal_max, "fix_nve_gpu:dtfm");
      }
    }
  }

  #if (LAL_USE_OMP == 1)
  #pragma omp parallel
  #endif
  {
    #if (LAL_USE_OMP == 1)
    const int nthreads = omp_get_num_threads();
    const int tid = omp_get_thread_num();
    const int idelta = nlocal / nthreads + 1;
    const int ifrom = tid * idelta;
    const int ito = MIN(ifrom + idelta, nlocal);
    const int ifrom3 = ifrom * 3;
    const int ito3 = ito * 3;
    #else
    const int tid = 0;
    const int ifrom = 0;
    const int ifrom3 = 0;
    const int ito = nlocal;
    const int ito3 = nlocal * 3;
    #endif
    if (igroup == 0 && atom->ntypes == 1 && !atom->rmass) {
      const double dtfm = dtf / atom->mass[1];
      #if (LAL_USE_OMP_SIMD == 1)
      #pragma omp simd
      #endif
      for (int i = ifrom3; i < ito3; i++)
        v[i] += dtfm * f[i];
    } else {
      if (neighbor->ago == 0) reset_dt_omp(ifrom,ito,tid);
      for (int i = ifrom; i < ito; i++) {
        const double dtfm = _dtfm[i];
        const int n = i * 3;
        v[n] += dtfm * f[n];
        v[n+1] += dtfm * f[n+1];
        v[n+2] += dtfm * f[n+2];
      }
    }
  }
}

void FixNVEGPU::reset_dt() {
  if (_respa_on) { FixNVE::reset_dt(); return; }
  if (igroup == 0 && atom->ntypes == 1 && !atom->rmass) {
    dtv = update->dt;
    dtf = 0.5 * update->dt * force->ftm2v;
  } else {
    const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
      atom->nlocal;
    if (nlocal > _nlocal_max) {
      if (_nlocal_max) memory->destroy(_dtfm);
      _nlocal_max = static_cast<int>(1.20 * nlocal);
      memory->create(_dtfm, _nlocal_max, "fix_nve_gpu:dtfm");
    }

    #if (LAL_USE_OMP == 1)
    #pragma omp parallel
    #endif
    {
      #if (LAL_USE_OMP == 1)
      const int nthreads = omp_get_num_threads();
      const int tid = omp_get_thread_num();
      const int idelta = nlocal / nthreads + 1;
      const int ifrom = tid * idelta;
      const int ito = MIN(ifrom + idelta, nlocal);
      #else
      const int tid = 0;
      const int ifrom = 0;
      const int ito = nlocal;
      #endif

      reset_dt_omp(ifrom, ito, tid);
    }
  }
}

void FixNVEGPU::reset_dt_omp(const int ifrom, const int ito, const int tid) {
  const double dtfo = 0.5 * update->dt * force->ftm2v;
  if (tid == 0) {
    dtv = update->dt;
    dtf = dtfo;
  }

  const int * const mask = atom->mask;
  if (igroup == 0) {
    if (atom->rmass) {
      const double * const rmass = atom->rmass;
      for (int i = ifrom; i < ito; i++) _dtfm[i] = dtfo / rmass[i];
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      for (int i = ifrom; i < ito; i++) _dtfm[i] = dtfo / mass[type[i]];
    }
  } else {
    if (atom->rmass) {
      const double * const rmass = atom->rmass;
      for (int i = ifrom; i < ito; i++)
        _dtfm[i] = (mask[i] & groupbit) ? dtfo / rmass[i] : 0.0;
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      for (int i = ifrom; i < ito; i++)
        _dtfm[i] = (mask[i] & groupbit) ? dtfo / mass[type[i]] : 0.0;
    }
  }
}

double FixNVEGPU::memory_usage()
{
  return FixNVE::memory_usage() + _nlocal_max * sizeof(double);
}
