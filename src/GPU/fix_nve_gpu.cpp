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
#include "comm.h"
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
  const int nlocal3 = nlocal * 3;

  #if (LAL_USE_OMP == 1)
  #pragma omp parallel
  #endif
  {
    #if (LAL_USE_OMP == 1)
    const int nthreads = comm->nthreads;
    const int idelta = nlocal3 / nthreads + 1;
    const int ifrom3 = omp_get_thread_num() * idelta;
    const int ito3 = MIN(ifrom3 + idelta, nlocal3);
    #else
    const int ifrom3 = 0;
    const int ito3 = nlocal3;
    #endif
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
      #if (LAL_USE_OMP_SIMD == 1)
      #pragma omp simd
      #endif
      for (int i = ifrom3; i < ito3; i++) {
        v[i] += _dtfm[i] * f[i];
        x[i] += dtv * v[i];
      }
    } else {
      #if (LAL_USE_OMP_SIMD == 1)
      #pragma omp simd
      #endif
      for (int i = ifrom3; i < ito3; i++) {
        if (_dtfm[i] != 0.0) {
          v[i] += _dtfm[i] * f[i];
          x[i] += dtv * v[i];
        }
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
        memory->create(_dtfm, _nlocal_max * 3, "fix_nve_gpu:dtfm");
      }
    }
  }

  #if (LAL_USE_OMP == 1)
  #pragma omp parallel
  #endif
  {
    #if (LAL_USE_OMP == 1)
    const int nthreads = comm->nthreads;
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
    } else if (igroup == 0) {
      if (neighbor->ago == 0) reset_dt_omp(ifrom,ito,tid);
      #if (LAL_USE_OMP_SIMD == 1)
      #pragma omp simd
      #endif
      for (int i = ifrom3; i < ito3; i++)
        v[i] += _dtfm[i] * f[i];
    } else {
      if (neighbor->ago == 0) reset_dt_omp(ifrom,ito,tid);
      #if (LAL_USE_OMP_SIMD == 1)
      #pragma omp simd
      #endif
      for (int i = ifrom3; i < ito3; i++)
        v[i] += _dtfm[i] * f[i];
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
      memory->create(_dtfm, _nlocal_max * 3, "fix_nve_gpu:dtfm");
    }

    #if (LAL_USE_OMP == 1)
    #pragma omp parallel
    #endif
    {
      #if (LAL_USE_OMP == 1)
      const int nthreads = comm->nthreads;
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
      int n = ifrom * 3;
      for (int i = ifrom; i < ito; i++) {
        const double dtfir = dtfo / rmass[i];
        _dtfm[n++] = dtfir;
        _dtfm[n++] = dtfir;
        _dtfm[n++] = dtfir;
      }
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      int n = ifrom * 3;
      for (int i = ifrom; i < ito; i++) {
        const double dtfim = dtfo / mass[type[i]];
        _dtfm[n++] = dtfim;
        _dtfm[n++] = dtfim;
        _dtfm[n++] = dtfim;
      }
    }
  } else {
    if (atom->rmass) {
      const double * const rmass = atom->rmass;
      int n = ifrom * 3;
      for (int i = ifrom; i < ito; i++)
        if (mask[i] & groupbit) {
          const double dtfir = dtfo / rmass[i];
          _dtfm[n++] = dtfir;
          _dtfm[n++] = dtfir;
          _dtfm[n++] = dtfir;
        } else {
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
        }
    } else {
      const double * const mass = atom->mass;
      const int * const type = atom->type;
      int n = ifrom * 3;
      for (int i = ifrom; i < ito; i++)
        if (mask[i] & groupbit) {
          const double dtfim = dtfo / mass[type[i]];
          _dtfm[n++] = dtfim;
          _dtfm[n++] = dtfim;
          _dtfm[n++] = dtfim;
        } else {
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
          _dtfm[n++] = 0.0;
        }
    }
  }
}

double FixNVEGPU::memory_usage()
{
  const int nlocal = (igroup == atom->firstgroup) ? atom->nfirst :
    atom->nlocal;
  return FixNVE::memory_usage() + nlocal * 3 * sizeof(double);
}
