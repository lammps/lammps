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

#include "wire_dipole_intel.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"

using namespace LAMMPS_NS;
using namespace MathConst;

WireDipoleIntel::WireDipoleIntel(LAMMPS *lmp, FixIntel *fix) : BoundaryCorrection(lmp)
{
  this->fix = fix;
  _use_lrt = fix->lrt();
}

/* ----------------------------------------------------------------------
   Wire-geometry correction term to dampen inter-wire interactions between
   periodically repeating wires.  Yields good approximation to 1D Ewald if
   adequate empty space is left between repeating wires (J. Mol. Struct.
   704, 101). x and y are non-periodic.
------------------------------------------------------------------------- */

void WireDipoleIntel::compute_corr(double qsum, int eflag_atom, int eflag_global, double &energy,
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
void WireDipoleIntel::compute_corr(IntelBuffers<flt_t, acc_t> *buffers, double /*qsum*/,
                                   int eflag_atom, int eflag_global, double &energy, double *eatom)
{
  double const volume = get_volume();
  ATOM_T *_noalias const x = buffers->get_x(0);
  flt_t *_noalias const q = buffers->get_q(0);
  int nlocal = atom->nlocal;
  int nthr;
  if (_use_lrt)
    nthr = 1;
  else
    nthr = comm->nthreads;

  double xdipole = 0.0;
  double ydipole = 0.0;
#if defined(_OPENMP)
#pragma omp parallel reduction(+ : xdipole, ydipole) \
    shared(nlocal, nthr, x, q) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      xdipole += q[i] * x[i].x;
      ydipole += q[i] * x[i].y;
    }
  }

  double xdipole_all, ydipole_all;
  MPI_Allreduce(&xdipole, &xdipole_all, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&ydipole, &ydipole_all, 1, MPI_DOUBLE, MPI_SUM, world);

  double xdipole_r2 = 0.0;
  double ydipole_r2 = 0.0;
  if (eflag_atom) {
#if defined(_OPENMP)
#pragma omp parallel reduction(+ : xdipole_r2, ydipole_r2) \
    shared(nlocal, nthr, x, q) if (!_use_lrt)
#endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

      for (int i = ifrom; i < ito; i++) {
        xdipole_r2 += q[i] * x[i].x * x[i].x;
        ydipole_r2 += q[i] * x[i].y * x[i].y;
      }
    }

    double tmp;
    MPI_Allreduce(&xdipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
    xdipole_r2 = tmp;
    MPI_Allreduce(&ydipole_r2, &tmp, 1, MPI_DOUBLE, MPI_SUM, world);
    ydipole_r2 = tmp;
  }

  const double e_wirecorr =
      MY_PI * (xdipole_all * xdipole_all + ydipole_all * ydipole_all) / volume;
  double const scale = 1.0;
  const double qscale = force->qqrd2e * scale;
  if (eflag_global) energy += qscale * e_wirecorr;

  if (eflag_atom) {
    double efact = qscale * MY_PI / volume;
#if defined(_OPENMP)
#pragma omp parallel shared(nlocal, nthr, x, q, eatom, efact, xdipole_all, \
                                                 ydipole_all, xdipole_r2,                   \
                                                 ydipole_r2) if (!_use_lrt)
#endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

      for (int i = ifrom; i < ito; i++)
        eatom[i] += efact * q[i] *
            (x[i].x * xdipole_all + x[i].y * ydipole_all - 0.5 * (xdipole_r2 + ydipole_r2));
    }
  }

  double ffact = qscale * (-MY_2PI / volume);
  double **f = atom->f;
#if defined(_OPENMP)
#pragma omp parallel shared(nlocal, nthr, x, q, f, ffact, xdipole_all, \
                                                 ydipole_all) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      f[i][0] += ffact * q[i] * xdipole_all;
      f[i][1] += ffact * q[i] * ydipole_all;
    }
  }
}

/* ----------------------------------------------------------------------
   vector_corr -- templates and implementation
------------------------------------------------------------------------- */

void WireDipoleIntel::vector_corr(double *vec, int sensor_grpbit, int source_grpbit,
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
void WireDipoleIntel::vector_corr(IntelBuffers<flt_t, acc_t> *buffers, double *vec,
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

  double xdipole = 0.0;
  double ydipole = 0.0;
#if defined(_OPENMP)
#pragma omp parallel reduction(+ : xdipole, ydipole) \
    shared(nlocal, nthr, x, q, mask, source_grpbit, invert_source) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (!!(mask[i] & source_grpbit) != invert_source) {
        xdipole += q[i] * x[i].x;
        ydipole += q[i] * x[i].y;
      }
    }
  }

  double dipole[2] = {xdipole, ydipole};
  MPI_Allreduce(MPI_IN_PLACE, &dipole, 2, MPI_DOUBLE, MPI_SUM, world);
  for (double &d : dipole) d *= MY_2PI / volume;

#if defined(_OPENMP)
#pragma omp parallel shared(nlocal, nthr, vec, x, mask, sensor_grpbit, \
                                                 dipole) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (mask[i] & sensor_grpbit) vec[i] += x[i].x * dipole[0] + x[i].y * dipole[1];
    }
  }
}

/* ----------------------------------------------------------------------
   matrix_corr -- templates and implementation
------------------------------------------------------------------------- */

void WireDipoleIntel::matrix_corr(bigint *imat, double **matrix)
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
void WireDipoleIntel::matrix_corr(IntelBuffers<flt_t, acc_t> *buffers, bigint *imat,
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
#pragma omp parallel reduction(+ : ngrouplocal) \
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

  std::vector<double> xprd_local = std::vector<double>(ngrouplocal);
  std::vector<double> yprd_local = std::vector<double>(ngrouplocal);

  std::vector<int> ngrouplocal_per_thr = std::vector<int>(nthr, 0);
#if defined(_OPENMP)
#pragma omp parallel shared(nlocal, nthr, imat, x, xprd_local, yprd_local, \
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
#pragma omp parallel shared(nlocal, nthr, imat, x, xprd_local, yprd_local, \
                                                 ngrouplocal_per_thr) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    int n = ngrouplocal_per_thr[tid];
    for (int i = ifrom; i < ito; i++) {
      if (imat[i] < 0) continue;
      xprd_local[n] = x[i].x;
      yprd_local[n] = x[i].y;
      n++;
    }
  }

  std::vector<int> recvcounts = gather_recvcounts(ngrouplocal);
  std::vector<int> displs = gather_displs(recvcounts);
  std::vector<double> xprd_all = std::vector<double>(ngroup);
  std::vector<double> yprd_all = std::vector<double>(ngroup);
  MPI_Allgatherv(xprd_local.data(), ngrouplocal, MPI_DOUBLE, xprd_all.data(), recvcounts.data(),
                 displs.data(), MPI_DOUBLE, world);
  MPI_Allgatherv(yprd_local.data(), ngrouplocal, MPI_DOUBLE, yprd_all.data(), recvcounts.data(),
                 displs.data(), MPI_DOUBLE, world);

  std::vector<bigint> jmat = gather_jmat(imat);
  const double prefac = MY_2PI / volume;
#if defined(_OPENMP)
#pragma omp parallel shared(nlocal, nthr, imat, x, jmat, ngroup, xprd_all, \
                                                 yprd_all, matrix, prefac) if (!_use_lrt)
#endif
  {
    int ifrom, ito, tid;
    IP_PRE_omp_range_id(ifrom, ito, tid, nlocal, nthr);

    for (int i = ifrom; i < ito; i++) {
      if (imat[i] < 0) continue;
      for (bigint j = 0; j < ngroup; j++) {
        if (jmat[j] > imat[i]) continue;
        double aij = prefac * (x[i].x * xprd_all[j] + x[i].y * yprd_all[j]);
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
