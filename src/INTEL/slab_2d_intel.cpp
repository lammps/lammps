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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Robert Meissner (Hereon, TUHH),
   Shern Tee (GU) with LLM (GLM-5.1)
------------------------------------------------------------------------- */

#include "slab_2d_intel.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

Slab2dIntel::Slab2dIntel(LAMMPS *lmp, FixIntel *fix) : BoundaryCorrection(lmp)
{
  this->fix = fix;
  _use_lrt = fix->lrt();
}

/* ----------------------------------------------------------------------
   Slab-geometry correction term (k=0) of EW2D. See Hu, JCTC 10:12 (2014)
   pp. 5254-5264 or metalwalls ewald and parallelization documentation.
------------------------------------------------------------------------- */

void Slab2dIntel::compute_corr(double qsum, int eflag_atom, int eflag_global, double &energy,
                               double *eatom)
{
  switch (fix->precision()) {
    case FixIntel::PREC_MODE_MIXED:
      compute_corr<float, double>(fix->get_mixed_buffers(), qsum, eflag_atom, eflag_global, energy,
                                  eatom);
      break;
    case FixIntel::PREC_MODE_DOUBLE:
      compute_corr<double, double>(fix->get_double_buffers(), qsum, eflag_atom, eflag_global,
                                   energy, eatom);
      break;
    default:
      compute_corr<float, float>(fix->get_single_buffers(), qsum, eflag_atom, eflag_global, energy,
                                 eatom);
  }
}

template <class flt_t, class acc_t>
void Slab2dIntel::compute_corr(IntelBuffers<flt_t, acc_t> *buffers, double /*qsum*/, int eflag_atom,
                               int eflag_global, double &energy, double *eatom)
{
  flt_t *_noalias const q = buffers->get_q();
  ATOM_T *_noalias const x = buffers->get_x();
  double **f = atom->f;
  int nlocal = atom->nlocal;
  double const g_ewald = force->kspace->g_ewald;
  bigint natoms = atom->natoms;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  std::vector<double> z = std::vector<double>(nlocal);

  #if defined(_OPENMP)
  #pragma omp parallel \
    shared(nlocal, nthr, x, z) if (!_use_lrt)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) z[i] = x[i].z;
  }

  std::vector<double> z_all = std::vector<double>(natoms);
  std::vector<double> q_all = std::vector<double>(natoms);
  std::vector<int> recvcounts = gather_recvcounts(nlocal);
  std::vector<int> displs = gather_displs(recvcounts);
  MPI_Allgatherv(q, nlocal, MPI_DOUBLE, q_all.data(), recvcounts.data(), displs.data(), MPI_DOUBLE,
                 world);
  MPI_Allgatherv(z.data(), nlocal, MPI_DOUBLE, z_all.data(), recvcounts.data(), displs.data(),
                 MPI_DOUBLE, world);

  const double g_ewald_inv = 1.0 / g_ewald;
  double const scale = 1.0;
  const double qscale = force->qqrd2e * scale;
  double const area = domain->xprd * domain->yprd;
  const double ffact = qscale * MY_2PI / area;
  const double efact = qscale * MY_PIS / area;

  double e_keq0 = 0;

  #if defined(_OPENMP)
  #pragma omp parallel \
    shared(nlocal, nthr, x, q, f, q_all, z_all, natoms, g_ewald, g_ewald_inv, ffact, efact, \
    eflag_atom, eflag_global, eatom) \
    reduction(+ : e_keq0) if (!_use_lrt)
  #endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      double pot_ij = 0.0;

      for (bigint j = 0; j < natoms; j++) {
        double const zij = z_all[j] - x[i].z;
        double const g_zij = g_ewald * zij;

        pot_ij += q_all[j] * (exp(-g_zij * g_zij) * g_ewald_inv + MY_PIS * zij * erf(g_zij));
        f[i][2] -= ffact * q[i] * q_all[j] * erf(g_zij);
      }

      if (eflag_atom) eatom[i] -= efact * q[i] * pot_ij;
      if (eflag_global) e_keq0 -= q[i] * pot_ij;
    }
  }

  if (eflag_global) {
    MPI_Allreduce(MPI_IN_PLACE, &e_keq0, 1, MPI_DOUBLE, MPI_SUM, world);
    energy += efact * e_keq0;
  }
}

/* ----------------------------------------------------------------------
   vector_corr -- templates and implementation
------------------------------------------------------------------------- */

void Slab2dIntel::vector_corr(double *vec, int sensor_grpbit, int source_grpbit, bool invert_source)
{
  switch (fix->precision()) {
    case FixIntel::PREC_MODE_MIXED:
      vector_corr<float, double>(fix->get_mixed_buffers(), vec, sensor_grpbit, source_grpbit,
                                 invert_source);
      break;
    case FixIntel::PREC_MODE_DOUBLE:
      vector_corr<double, double>(fix->get_double_buffers(), vec, sensor_grpbit, source_grpbit,
                                  invert_source);
      break;
    default:
      vector_corr<float, float>(fix->get_single_buffers(), vec, sensor_grpbit, source_grpbit,
                                invert_source);
  }
}

template <class flt_t, class acc_t>
void Slab2dIntel::vector_corr(IntelBuffers<flt_t, acc_t> *buffers, double *vec, int sensor_grpbit,
                              int source_grpbit, bool invert_source)
{
  int const nlocal = atom->nlocal;
  ATOM_T *_noalias const x = buffers->get_x();
  flt_t *_noalias const q = buffers->get_q();
  int *mask = atom->mask;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  std::vector<double> z_local;
  std::vector<double> q_local;

  int nelectrolyte_local = 0;
#if defined(_OPENMP)
#pragma omp parallel \
 reduction(+ : nelectrolyte_local) \
    shared(nlocal, nthr, mask, source_grpbit, invert_source) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++)
      if (!!(mask[i] & source_grpbit) != invert_source) nelectrolyte_local++;
  }

  z_local.resize(nelectrolyte_local);
  q_local.resize(nelectrolyte_local);

  std::vector<int> nelectrolyte_per_thr = std::vector<int>(nthr, 0);
#if defined(_OPENMP)
#pragma omp parallel \
 shared(nlocal, nthr, mask, source_grpbit, invert_source, x, \
                                                 q, z_local, q_local,                             \
                                                 nelectrolyte_per_thr) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    int n = 0;
    for (int i = ifrom; i < ito; i++)
      if (!!(mask[i] & source_grpbit) != invert_source) n++;
    nelectrolyte_per_thr[tid] = n;
  }

  int offset = 0;
  for (int t = 0; t < nthr; t++) {
    int tmp = nelectrolyte_per_thr[t];
    nelectrolyte_per_thr[t] = offset;
    offset += tmp;
  }

#if defined(_OPENMP)
#pragma omp parallel \
 shared(nlocal, nthr, mask, source_grpbit, invert_source, x, \
                                                 q, z_local, q_local,                             \
                                                 nelectrolyte_per_thr) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    int n = nelectrolyte_per_thr[tid];
    for (int i = ifrom; i < ito; i++) {
      if (!!(mask[i] & source_grpbit) != invert_source) {
        z_local[n] = x[i].z;
        q_local[n] = q[i];
        n++;
      }
    }
  }

  int n_electrolyte;
  MPI_Allreduce(&nelectrolyte_local, &n_electrolyte, 1, MPI_INT, MPI_SUM, world);
  std::vector<double> z_all = std::vector<double>(n_electrolyte);
  std::vector<double> q_all = std::vector<double>(n_electrolyte);
  std::vector<int> recvcounts = gather_recvcounts(nelectrolyte_local);
  std::vector<int> displs = gather_displs(recvcounts);
  MPI_Allgatherv(z_local.data(), nelectrolyte_local, MPI_DOUBLE, z_all.data(), recvcounts.data(),
                 displs.data(), MPI_DOUBLE, world);
  MPI_Allgatherv(q_local.data(), nelectrolyte_local, MPI_DOUBLE, q_all.data(), recvcounts.data(),
                 displs.data(), MPI_DOUBLE, world);

  double const g_ewald = force->kspace->g_ewald;
  double const g_ewald_inv = 1.0 / g_ewald;
  double const area = domain->xprd * domain->yprd;
  double const prefac = 2 * MY_PIS / area;
#if defined(_OPENMP)
#pragma omp parallel \
 shared(nlocal, nthr, vec, x, mask, sensor_grpbit, z_all, \
                                                 q_all, n_electrolyte, g_ewald, g_ewald_inv,   \
                                                 prefac) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (!(mask[i] & sensor_grpbit)) continue;
      double b = 0;
      double zi = x[i].z;
      for (int j = 0; j < n_electrolyte; j++) {
        double zij = z_all[j] - zi;
        double gzij = g_ewald * zij;
        double zfac = (gzij > 7) ? MY_PIS * zij
                                 : exp(-(gzij * gzij)) / g_ewald + MY_PIS * zij * erf(gzij);
        b += q_all[j] * zfac;
      }
      vec[i] -= prefac * b;
    }
  }
}

/* ----------------------------------------------------------------------
   matrix_corr -- templates and implementation
------------------------------------------------------------------------- */

void Slab2dIntel::matrix_corr(bigint *imat, double **matrix)
{
  switch (fix->precision()) {
    case FixIntel::PREC_MODE_MIXED:
      matrix_corr<float, double>(fix->get_mixed_buffers(), imat, matrix);
      break;
    case FixIntel::PREC_MODE_DOUBLE:
      matrix_corr<double, double>(fix->get_double_buffers(), imat, matrix);
      break;
    default:
      matrix_corr<float, float>(fix->get_single_buffers(), imat, matrix);
  }
}

template <class flt_t, class acc_t>
void Slab2dIntel::matrix_corr(IntelBuffers<flt_t, acc_t> *buffers, bigint *imat, double **matrix)
{
  int nlocal = atom->nlocal;
  ATOM_T *_noalias const x = buffers->get_x();
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  int ngrouplocal = 0;
#if defined(_OPENMP)
#pragma omp parallel \
 reduction(+ : ngrouplocal) \
    shared(nlocal, nthr, imat) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++)
      if (imat[i] > -1) ngrouplocal++;
  }

  bigint ngroup = 0;
  MPI_Allreduce(&ngrouplocal, &ngroup, 1, MPI_INT, MPI_SUM, world);

  std::vector<double> nprd_local = std::vector<double>(ngrouplocal);
  std::vector<int> ngrouplocal_per_thr = std::vector<int>(nthr, 0);
#if defined(_OPENMP)
#pragma omp parallel \
 shared(nlocal, nthr, imat, x, nprd_local, \
                                                 ngrouplocal_per_thr) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    int n = 0;
    for (int i = ifrom; i < ito; i++)
      if (imat[i] > -1) n++;
    ngrouplocal_per_thr[tid] = n;
  }

  int offset = 0;
  for (int t = 0; t < nthr; t++) {
    int tmp = ngrouplocal_per_thr[t];
    ngrouplocal_per_thr[t] = offset;
    offset += tmp;
  }

#if defined(_OPENMP)
#pragma omp parallel \
 shared(nlocal, nthr, imat, x, nprd_local, \
                                                 ngrouplocal_per_thr) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    int n = ngrouplocal_per_thr[tid];
    for (int i = ifrom; i < ito; i++) {
      if (imat[i] < 0) continue;
      nprd_local[n++] = x[i].z;
    }
  }

  std::vector<int> recvcounts = gather_recvcounts(ngrouplocal);
  std::vector<int> displs = gather_displs(recvcounts);
  std::vector<double> nprd_all = std::vector<double>(ngroup);
  MPI_Allgatherv(nprd_local.data(), ngrouplocal, MPI_DOUBLE, nprd_all.data(), recvcounts.data(),
                 displs.data(), MPI_DOUBLE, world);

  double const g_ewald = force->kspace->g_ewald;
  const double g_ewald_inv = 1.0 / g_ewald;
  const double g_ewald_sq = g_ewald * g_ewald;
  double const area = domain->xprd * domain->yprd;
  const double prefac = 2.0 * MY_PIS / area;
  std::vector<bigint> jmat = gather_jmat(imat);
#if defined(_OPENMP)
#pragma omp parallel \
 shared(nlocal, nthr, imat, x, jmat, ngroup, nprd_all, \
                                                 matrix, prefac, g_ewald, g_ewald_inv,      \
                                                 g_ewald_sq) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (imat[i] < 0) continue;
      for (bigint j = 0; j < ngroup; j++) {
        if (jmat[j] > imat[i]) continue;
        double dij = nprd_all[j] - x[i].z;
        double aij = prefac *
            (exp(-dij * dij * g_ewald_sq) * g_ewald_inv + MY_PIS * dij * erf(dij * g_ewald));
#if defined(_OPENMP)
#pragma omp atomic
#endif
        matrix[imat[i]][jmat[j]] -= aij;
        if (imat[i] != jmat[j]) {
#if defined(_OPENMP)
#pragma omp atomic
#endif
          matrix[jmat[j]][imat[i]] -= aij;
        }
      }
    }
  }
}
