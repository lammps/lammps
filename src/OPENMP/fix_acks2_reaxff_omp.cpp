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

#include "fix_acks2_reaxff_omp.h"

#include "atom.h"
#include "comm.h"
#include "memory.h"

#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixACKS2ReaxFFOMP::FixACKS2ReaxFFOMP(LAMMPS *lmp, int narg, char **arg) :
  FixACKS2ReaxFF(lmp, narg, arg), b_temp(nullptr), nmax_temp(0)
{
}

/* ---------------------------------------------------------------------- */

FixACKS2ReaxFFOMP::~FixACKS2ReaxFFOMP()
{
  memory->destroy(b_temp);
}

/* ---------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::allocate_storage()
{
  FixACKS2ReaxFF::allocate_storage();

  // per-thread buffers for the contributions to the H and X parts of the
  // matrix-vector product and to the diagonal of the X matrix

  nmax_temp = 2*nmax;
  memory->create(b_temp, comm->nthreads, nmax_temp, "acks2/reaxff/omp:b_temp");
}

/* ---------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::deallocate_storage()
{
  memory->destroy(b_temp);
  nmax_temp = 0;

  FixACKS2ReaxFF::deallocate_storage();
}

/* ---------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::compute_H()
{
  compute_H_offsets();

#if defined(_OPENMP)
#pragma omp parallel for schedule(guided)
#endif
  for (int ii = 0; ii < nn; ii++) compute_H_row(ii);
}

/* ----------------------------------------------------------------------
   the diagonal of the X matrix is accumulated in per-thread buffers that
   are summed in a fixed order of the threads
------------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::compute_X()
{
  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;

  compute_X_offsets();

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    double *xdiag = b_temp[tid];
    for (int i = 0; i < nall; ++i) xdiag[i] = 0.0;

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(guided)
#endif
    for (int ii = 0; ii < nn; ii++) compute_X_row(ii, xdiag);

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < nall; ++i) {
      double sum = 0.0;
      for (int t = 0; t < nthreads; ++t) sum += b_temp[t][i];
      X_diag[i] = sum;
    }
  }
}

/* ----------------------------------------------------------------------
   threaded version of FixACKS2ReaxFF::sparse_matvec_acks2().  each row is
   processed by one thread, the contributions to the other atom of a pair
   are collected in per-thread buffers and summed in a fixed order.
------------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::sparse_matvec_acks2(sparse_matrix *H, sparse_matrix *X, double *x,
                                            double *b)
{
  const int nlocal = atom->nlocal;
  const int nthreads = comm->nthreads;
  const int *mask = atom->mask;
  const int *type = atom->type;
  double last0 = 0.0, last1 = 0.0;

#if defined(_OPENMP)
#pragma omp parallel default(shared) reduction(+:last0,last1)
#endif
  {
#if defined(_OPENMP)
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    double *bt = b_temp[tid];
    for (int i = 0; i < 2*NN; ++i) bt[i] = 0.0;

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (int ii = 0; ii < nn; ++ii) {
      const int i = ilist[ii];
      if (mask[i] & groupbit) {
        b[i] = eta[type[i]] * x[i];
        b[NN + i] = X_diag[i] * x[NN + i];
      } else {
        b[i] = b[NN + i] = 0.0;
      }
    }

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (int i = nlocal; i < NN; ++i) {
      b[i] = 0.0;
      b[NN + i] = 0.0;
    }

    // the rows of atoms outside the group only contain neighbors in the
    // group for the H matrix and are empty for the X matrix

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (int ii = 0; ii < nn; ++ii) {
      const int i = ilist[ii];

      for (int itr_j = H->firstnbr[i]; itr_j < H->firstnbr[i] + H->numnbrs[i]; itr_j++) {
        const int j = H->jlist[itr_j];
        b[i] += H->val[itr_j] * x[j];
        bt[j] += H->val[itr_j] * x[i];
      }

      if (mask[i] & groupbit) {
        for (int itr_j = X->firstnbr[i]; itr_j < X->firstnbr[i] + X->numnbrs[i]; itr_j++) {
          const int j = X->jlist[itr_j];
          b[NN + i] += X->val[itr_j] * x[NN + j];
          bt[NN + j] += X->val[itr_j] * x[NN + i];
        }

        // identity matrix
        b[NN + i] += x[i];
        b[i] += x[NN + i];

        // second-to-last and last row/column
        last0 += x[NN + i];
        b[NN + i] += x[2*NN];
        last1 += x[i];
        b[i] += x[2*NN + 1];
      }
    }

#if defined(_OPENMP)
#pragma omp for schedule(static)
#endif
    for (int i = 0; i < 2*NN; ++i)
      for (int t = 0; t < nthreads; ++t) b[i] += b_temp[t][i];
  }

  b[2*NN] = last0;
  b[2*NN + 1] = last1;
}

/* ---------------------------------------------------------------------- */

double FixACKS2ReaxFFOMP::parallel_norm(double *v, int n)
{
  const int *mask = atom->mask;
  double my_sum = 0.0, norm_sqr = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) reduction(+:my_sum)
#endif
  for (int ii = 0; ii < n; ++ii) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      my_sum += v[i] * v[i];
      my_sum += v[NN + i] * v[NN + i];
    }
  }

  // last two rows
  if (last_rows_flag) {
    my_sum += v[2*NN] * v[2*NN];
    my_sum += v[2*NN + 1] * v[2*NN + 1];
  }

  MPI_Allreduce(&my_sum, &norm_sqr, 1, MPI_DOUBLE, MPI_SUM, world);

  return sqrt(norm_sqr);
}

/* ---------------------------------------------------------------------- */

double FixACKS2ReaxFFOMP::parallel_dot(double *v1, double *v2, int n)
{
  const int *mask = atom->mask;
  double my_dot = 0.0, res = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) reduction(+:my_dot)
#endif
  for (int ii = 0; ii < n; ++ii) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      my_dot += v1[i] * v2[i];
      my_dot += v1[NN + i] * v2[NN + i];
    }
  }

  // last two rows
  if (last_rows_flag) {
    my_dot += v1[2*NN] * v2[2*NN];
    my_dot += v1[2*NN + 1] * v2[2*NN + 1];
  }

  MPI_Allreduce(&my_dot, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

double FixACKS2ReaxFFOMP::parallel_vector_acc(double *v, int n)
{
  const int *mask = atom->mask;
  double my_acc = 0.0, res = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static) reduction(+:my_acc)
#endif
  for (int ii = 0; ii < n; ++ii) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      my_acc += v[i];
      my_acc += v[NN + i];
    }
  }

  // last two rows
  if (last_rows_flag) {
    my_acc += v[2*NN];
    my_acc += v[2*NN + 1];
  }

  MPI_Allreduce(&my_acc, &res, 1, MPI_DOUBLE, MPI_SUM, world);

  return res;
}

/* ---------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::vector_sum(double *dest, double c, double *v, double d, double *y, int k)
{
  const int *mask = atom->mask;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int ii = 0; ii < k; ++ii) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      dest[i] = c * v[i] + d * y[i];
      dest[NN + i] = c * v[NN + i] + d * y[NN + i];
    }
  }

  // last two rows
  if (last_rows_flag) {
    dest[2*NN] = c * v[2*NN] + d * y[2*NN];
    dest[2*NN + 1] = c * v[2*NN + 1] + d * y[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::vector_add(double *dest, double c, double *v, int k)
{
  const int *mask = atom->mask;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int ii = 0; ii < k; ++ii) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      dest[i] += c * v[i];
      dest[NN + i] += c * v[NN + i];
    }
  }

  // last two rows
  if (last_rows_flag) {
    dest[2*NN] += c * v[2*NN];
    dest[2*NN + 1] += c * v[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

void FixACKS2ReaxFFOMP::vector_copy(double *dest, double *v, int k)
{
  const int *mask = atom->mask;

#if defined(_OPENMP)
#pragma omp parallel for schedule(static)
#endif
  for (int ii = 0; ii < k; ++ii) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      dest[i] = v[i];
      dest[NN + i] = v[NN + i];
    }
  }

  // last two rows
  if (last_rows_flag) {
    dest[2*NN] = v[2*NN];
    dest[2*NN + 1] = v[2*NN + 1];
  }
}

/* ---------------------------------------------------------------------- */

double FixACKS2ReaxFFOMP::memory_usage()
{
  double bytes = FixACKS2ReaxFF::memory_usage();
  bytes += (double)comm->nthreads * nmax_temp * sizeof(double);
  return bytes;
}
