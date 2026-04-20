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

#include "slab_dipole_intel.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 0.00001;

SlabDipoleIntel::SlabDipoleIntel(LAMMPS *lmp, FixIntel *fix) : BoundaryCorrection(lmp)
{
  this->fix = fix;
  _use_lrt = fix->lrt();
};

/* ----------------------------------------------------------------------
   Slab-geometry correction term to dampen inter-slab interactions between
   periodically repeating slabs.  Yields good approximation to 2D Ewald if
   adequate empty space is left between repeating slabs (J. Chem. Phys.
   111, 3155).  Slabs defined here to be parallel to the xy plane. Also
   extended to non-neutral systems (J. Chem. Phys. 131, 094107).
-------------------------------------------------------------------------
*/

void SlabDipoleIntel::compute_corr(double qsum, int eflag_atom, int eflag_global, double &energy,
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
void SlabDipoleIntel::compute_corr(IntelBuffers<flt_t, acc_t> *buffers, double qsum, int eflag_atom,
                                   int eflag_global, double &energy, double *eatom)
{
  double const volume = get_volume();
  ATOM_T *_noalias const x = buffers->get_x(0);
  flt_t *_noalias const q = buffers->get_q(0);
  double zprd_slab = domain->zprd * force->kspace->slab_volfactor;
  int nlocal = atom->nlocal;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  double dipole = 0.0;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE reduction(+ : dipole) \
    shared(nlocal, nthr, x, q) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) dipole += q[i] * x[i].z;
  }

  MPI_Allreduce(MPI_IN_PLACE, &dipole, 1, MPI_DOUBLE, MPI_SUM, world);

  double dipole_r2 = 0.0;
  if (eflag_atom || fabs(qsum) > SMALL) {
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE reduction(+ : dipole_r2) \
    shared(nlocal, nthr, x, q) if (!_use_lrt)
#endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

      for (int i = ifrom; i < ito; i++) dipole_r2 += q[i] * x[i].z * x[i].z;
    }

    MPI_Allreduce(MPI_IN_PLACE, &dipole_r2, 1, MPI_DOUBLE, MPI_SUM, world);
  }

  double const e_slabcorr = MY_2PI *
      (dipole * dipole - qsum * dipole_r2 - qsum * qsum * zprd_slab * zprd_slab / 12.0) /
      volume;
  double const scale = 1.0;
  double const qscale = force->qqrd2e * scale;
  if (eflag_global) energy += qscale * e_slabcorr;

  if (eflag_atom) {
    double efact = qscale * MY_2PI / volume;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE shared(nlocal, nthr, x, q, eatom, efact, dipole, \
                                                 dipole_r2, qsum, zprd_slab) if (!_use_lrt)
#endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

      for (int i = ifrom; i < ito; i++)
        eatom[i] += efact * q[i] *
            (x[i].z * dipole - 0.5 * (dipole_r2 + qsum * x[i].z * x[i].z) -
             qsum * zprd_slab * zprd_slab / 12.0);
    }
  }

  double const ffact = qscale * (-4.0 * MY_PI / volume);
  double **f = atom->f;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE shared(nlocal, nthr, x, q, f, ffact, dipole, \
                                                 qsum) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) f[i][2] += ffact * q[i] * (dipole_all - qsum * x[i].z);
  }
}

/* ----------------------------------------------------------------------
   vector_corr -- templates and implementation
------------------------------------------------------------------------- */

void SlabDipoleIntel::vector_corr(double *vec, int sensor_grpbit, int source_grpbit,
                                  bool invert_source)
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
void SlabDipoleIntel::vector_corr(IntelBuffers<flt_t, acc_t> *buffers, double *vec,
                                  int sensor_grpbit, int source_grpbit, bool invert_source)
{
  double const volume = get_volume();
  int const nlocal = atom->nlocal;
  ATOM_T *_noalias const x = buffers->get_x(0);
  flt_t *_noalias const q = buffers->get_q(0);
  int *mask = atom->mask;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  double dipole = 0.0;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE reduction(+ : dipole) \
    shared(nlocal, nthr, x, q, mask, source_grpbit, invert_source) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (!!(mask[i] & source_grpbit) != invert_source) dipole += q[i] * x[i].z;
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, &dipole, 1, MPI_DOUBLE, MPI_SUM, world);
  dipole *= 4.0 * MY_PI / volume;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE shared(nlocal, nthr, vec, x, mask, sensor_grpbit, \
                                                 dipole) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (mask[i] & sensor_grpbit) vec[i] += x[i].z * dipole;
    }
  }
}

/* ----------------------------------------------------------------------
   matrix_corr -- templates and implementation
------------------------------------------------------------------------- */

void SlabDipoleIntel::matrix_corr(bigint *imat, double **matrix)
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
void SlabDipoleIntel::matrix_corr(IntelBuffers<flt_t, acc_t> *buffers, bigint *imat,
                                  double **matrix)
{
  double const volume = get_volume();
  int nlocal = atom->nlocal;
  ATOM_T *_noalias const x = buffers->get_x(0);
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  int ngrouplocal = 0;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE reduction(+ : ngrouplocal) \
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
#pragma omp parallel LMP_DEFAULT_NONE shared(nlocal, nthr, imat, x, nprd_local, \
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
#pragma omp parallel LMP_DEFAULT_NONE shared(nlocal, nthr, imat, x, nprd_local, \
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

  std::vector<bigint> jmat = gather_jmat(imat);
  const double prefac = MY_4PI / volume;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE shared(nlocal, nthr, imat, x, jmat, ngroup, nprd_all, \
                                                 matrix, prefac) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (imat[i] < 0) continue;
      for (bigint j = 0; j < ngroup; j++) {
        if (jmat[j] > imat[i]) continue;
        double aij = prefac * x[i].z * nprd_all[j];
#if defined(_OPENMP)
#pragma omp atomic
#endif
        matrix[imat[i]][jmat[j]] += aij;
        if (imat[i] != jmat[j]) {
#if defined(_OPENMP)
#pragma omp atomic
#endif
          matrix[jmat[j]][imat[i]] += aij;
        }
      }
    }
  }
}
