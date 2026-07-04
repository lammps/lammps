// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_qeq_slater_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
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
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixQEqSlaterOMP::FixQEqSlaterOMP(LAMMPS *lmp, int narg, char **arg) :
  FixQEqSlater(lmp, narg, arg), b_temp(nullptr), nmax_btmp(0)
{
}

/* ---------------------------------------------------------------------- */

FixQEqSlaterOMP::~FixQEqSlaterOMP()
{
  memory->destroy(b_temp);
}

/* ---------------------------------------------------------------------- */

void FixQEqSlaterOMP::pre_force(int /*vflag*/)
{
  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;

  if (atom->nmax > nmax) reallocate_storage();

  if (nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  const int nthreads = comm->nthreads;
  if (atom->nmax > nmax_btmp) {
    memory->destroy(b_temp);
    nmax_btmp = atom->nmax;
    memory->create(b_temp, nthreads, nmax_btmp, "qeq/slater/omp:b_temp");
  }

  init_matvec_thr();
  matvecs = CG(b_s, s);
  matvecs += CG(b_t, t);
  matvecs /= 2;
  calculate_Q();

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqSlaterOMP::init_matvec_thr()
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
   threaded build of the slater H matrix.  Two-pass (count -> serial
   sum-scan offsets -> fill) gives contiguous, non-overlapping CSR blocks
   that reproduce FixQEqSlater::compute_H().  The per-i chizj accumulation
   stays local to each atom's block, so no cross-atom race occurs.
------------------------------------------------------------------------- */

void FixQEqSlaterOMP::compute_H_thr()
{
  int *ilist, *numneigh, **firstneigh;
  int *type = atom->type;
  double **x = atom->x;

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
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    int cnt = 0;
    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      const double dx = x[i][0] - x[j][0];
      const double dy = x[i][1] - x[j][1];
      const double dz = x[i][2] - x[j][2];
      if (dx*dx + dy*dy + dz*dz <= cutoff_sq) cnt++;
    }
    H.numnbrs[i] = cnt;
  }

  // serial sum-scan: assign contiguous CSR offsets in ilist order

  m_fill = 0;
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    H.firstnbr[i] = m_fill;
    m_fill += H.numnbrs[i];
  }

  // pass 2: fill each atom's contiguous block (race-free)

#if defined(_OPENMP)
#pragma omp parallel for schedule(dynamic,50) default(shared)
#endif
  for (int ii = 0; ii < inum; ii++) {
    const int i = ilist[ii];
    const int itype = type[i];
    const double zei = zeta[itype];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    int mfill = H.firstnbr[i];
    double zjtmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      const int j = jlist[jj] & NEIGHMASK;
      const int jtype = type[j];
      const double zej = zeta[jtype];
      const double zj = zcore[jtype];

      const double dx = x[i][0] - x[j][0];
      const double dy = x[i][1] - x[j][1];
      const double dz = x[i][2] - x[j][2];
      const double rsq = dx*dx + dy*dy + dz*dz;
      if (rsq > cutoff_sq) continue;

      const double r = sqrt(rsq);
      H.jlist[mfill] = j;
      if (vtype == 0)
        H.val[mfill] = calculate_H(zei, zej, zj, r, zjtmp);
      else if (vtype == 1)
        H.val[mfill] = calculate_H_wolf(zei, zej, zj, r, zjtmp);
      else if (vtype == 2)
        H.val[mfill] = calculate_H_dsf(zei, zej, zj, r, zjtmp);
      mfill++;
    }
    chizj[i] = zjtmp;
  }

  if (m_fill >= H.m)
    error->all(FLERR,"Fix qeq/slater/omp has insufficient H "
                     "matrix size: m_fill={} H.m={}\n",m_fill,H.m);
}

/* ----------------------------------------------------------------------
   threaded sparse matrix-vector product (slater woself diagonal)
------------------------------------------------------------------------- */

void FixQEqSlaterOMP::sparse_matvec(sparse_matrix *A, double *x, double *b)
{
  const int nthreads = comm->nthreads;
  const double rc = cutoff;
  const double woself = 0.50*erfc(alpha*rc)/rc + alpha/MY_PIS;
  const double diag_shift = 2.0*force->qqr2e*woself;

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
  {
    int tid = 0;
#if defined(_OPENMP)
    tid = omp_get_thread_num();
#endif
    const int nloc = atom->nlocal;
    const int na = atom->nlocal + atom->nghost;
    double *b_t = b_temp[tid];

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (int i = 0; i < nloc; ++i)
      if (atom->mask[i] & groupbit) b[i] = (eta[atom->type[i]] - diag_shift) * x[i];

#if defined(_OPENMP)
#pragma omp for schedule(dynamic,50)
#endif
    for (int i = nloc; i < na; ++i)
      if (atom->mask[i] & groupbit) b[i] = 0.0;

    for (int i = 0; i < na; ++i) b_t[i] = 0.0;

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
    for (int i = 0; i < na; ++i)
      for (int t = 0; t < nthreads; ++t) b[i] += b_temp[t][i];
  }
}
