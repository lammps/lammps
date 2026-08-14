// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_qeq_point_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixQEqPointOMP::FixQEqPointOMP(LAMMPS *lmp, int narg, char **arg) :
  FixQEqPoint(lmp, narg, arg), b_temp(nullptr), nmax_btmp(0)
{
}

/* ---------------------------------------------------------------------- */

FixQEqPointOMP::~FixQEqPointOMP()
{
  memory->destroy(b_temp);
}

/* ---------------------------------------------------------------------- */

void FixQEqPointOMP::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;

  if (atom->nmax > nmax) reallocate_storage();

  if (nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  // (re)allocate per-thread scatter buffers for sparse_matvec

  const int nthreads = comm->nthreads;
  if (atom->nmax > nmax_btmp) {
    memory->destroy(b_temp);
    nmax_btmp = atom->nmax;
    memory->create(b_temp, nthreads, nmax_btmp, "qeq/point/omp:b_temp");
  }

  init_matvec_thr();
  matvecs = CG(b_s, s);         // CG on s - parallel
  matvecs += CG(b_t, t);        // CG on t - parallel
  matvecs /= 2;
  calculate_Q();

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqPointOMP::init_matvec_thr()
{
  compute_H_thr();

  int inum, ii, i;
  int *ilist;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {
      Hdia_inv[i] = 1. / eta[atom->type[i]];
      b_s[i]      = -(chi[atom->type[i]] + chizj[i]);
      b_t[i]      = -1.0;
      t[i] = t_hist[i][2] + 3 * (t_hist[i][0] - t_hist[i][1]);
      s[i] = 4*(s_hist[i][0]+s_hist[i][2])-(6*s_hist[i][1]+s_hist[i][3]);
    }
  }

  pack_flag = 2;
  comm->forward_comm(this); //Dist_vector(s);
  pack_flag = 3;
  comm->forward_comm(this); //Dist_vector(t);
}

/* ----------------------------------------------------------------------
   threaded build of the sparse H matrix.  A serial sum-scan over the
   per-atom (filtered) neighbor counts assigns each atom a contiguous,
   non-overlapping block in the CSR arrays, so the two per-atom fill passes
   are race-free.  The resulting matrix layout is identical to the serial
   FixQEqPoint::compute_H().
------------------------------------------------------------------------- */

void FixQEqPointOMP::compute_H_thr()
{
  int *ilist, *numneigh, **firstneigh;
  double **x = atom->x;
  int *mask = atom->mask;

  const int inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // pass 1: count neighbors within cutoff for each atom

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) default(shared)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      const int *jlist = firstneigh[i];
      const int jnum = numneigh[i];
      int cnt = 0;
      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const double dx = x[j][0] - x[i][0];
        const double dy = x[j][1] - x[i][1];
        const double dz = x[j][2] - x[i][2];
        if (dx*dx + dy*dy + dz*dz <= cutoff_sq) cnt++;
      }
      H.numnbrs[i] = cnt;
    }
  }

  // serial sum-scan: assign contiguous CSR offsets in ilist order

  m_fill = 0;
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      H.firstnbr[i] = m_fill;
      m_fill += H.numnbrs[i];
    }
  }

  // pass 2: fill each atom's contiguous block (race-free)

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) default(shared)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    if (mask[i] & groupbit) {
      const int *jlist = firstneigh[i];
      const int jnum = numneigh[i];
      int mfill = H.firstnbr[i];
      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const double dx = x[j][0] - x[i][0];
        const double dy = x[j][1] - x[i][1];
        const double dz = x[j][2] - x[i][2];
        const double r_sqr = dx*dx + dy*dy + dz*dz;
        if (r_sqr <= cutoff_sq) {
          H.jlist[mfill] = j;
          H.val[mfill] = 0.5/sqrt(r_sqr);
          mfill++;
        }
      }
    }
  }
}

/* ----------------------------------------------------------------------
   threaded sparse matrix-vector product.  The diagonal/gather part writes
   only b[i] (race-free under the omp for over i).  The symmetric scatter
   b[j] += val*x[i] is accumulated into per-thread buffers and reduced, to
   avoid races on shared j (ghosts / cross-thread neighbors).
------------------------------------------------------------------------- */

void FixQEqPointOMP::sparse_matvec(sparse_matrix *A, double *x, double *b)
{
  const int nthreads = comm->nthreads;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int tid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif
    const int nloc = atom->nlocal;
    const int nall = atom->nlocal + atom->nghost;
    double *b_t = b_temp[tid];

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (int i = 0; i < nloc; ++i)
      if (atom->mask[i] & groupbit) b[i] = eta[atom->type[i]] * x[i];

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (int i = nloc; i < nall; ++i)
      if (atom->mask[i] & groupbit) b[i] = 0.0;

    for (int i = 0; i < nall; ++i) b_t[i] = 0.0;

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (int i = 0; i < nloc; ++i) {
      if (atom->mask[i] & groupbit) {
        for (int itr_j = A->firstnbr[i]; itr_j < A->firstnbr[i]+A->numnbrs[i]; itr_j++) {
          int j = A->jlist[itr_j];
          b[i] += A->val[itr_j] * x[j];
          b_t[j] += A->val[itr_j] * x[i];
        }
      }
    }

#if defined(_OPENMP)
#pragma omp barrier
#pragma omp for schedule(dynamic,50)
#endif
    for (int i = 0; i < nall; ++i)
      for (int t = 0; t < nthreads; ++t) b[i] += b_temp[t][i];
  }
}
