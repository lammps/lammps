// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_qeq_ctip_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "ewald_const.h"
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
using namespace EwaldConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

FixQEqCTIPOMP::FixQEqCTIPOMP(LAMMPS *lmp, int narg, char **arg) :
  FixQEqCTIP(lmp, narg, arg), b_temp(nullptr), nmax_btmp(0)
{
}

/* ---------------------------------------------------------------------- */

FixQEqCTIPOMP::~FixQEqCTIPOMP()
{
  memory->destroy(b_temp);
}

/* ---------------------------------------------------------------------- */

void FixQEqCTIPOMP::pre_force(int /*vflag*/)
{
  int i, n, nout;

  if (update->ntimestep % nevery) return;

  nlocal = atom->nlocal;

  if (atom->nmax > nmax) reallocate_storage();

  if (nlocal > n_cap*DANGER_ZONE || m_fill > m_cap*DANGER_ZONE)
    reallocate_matrix();

  const int nthreads = comm->nthreads;
  if (atom->nmax > nmax_btmp) {
    memory->destroy(b_temp);
    nmax_btmp = atom->nmax;
    memory->create(b_temp, nthreads, nmax_btmp, "qeq/ctip/omp:b_temp");
  }

  for (i = 1; i <= maxrepeat; i++) {
    init_matvec_thr();
    matvecs = CG(b_s, s);
    matvecs += CG(b_t, t);
    matvecs /= 2;
    n = calculate_check_Q();
    MPI_Allreduce(&n, &nout, 1, MPI_INT, MPI_SUM, world);
    if (nout == 0) break;
  }

  if (i > maxrepeat && comm->me == 0)
    error->all(FLERR, Error::NOLASTLINE, "Fix qeq some charges not bound within the domain");

  if (force->kspace) force->kspace->qsum_qsq();
}

/* ---------------------------------------------------------------------- */

void FixQEqCTIPOMP::init_matvec_thr()
{
  compute_H_thr();

  int inum, ii, i;
  int *ilist;
  double *q = atom->q, qi;
  int *type = atom->type;

  inum = list->inum;
  ilist = list->ilist;

  for (ii = 0; ii < inum; ++ii) {
    i = ilist[ii];
    if (atom->mask[i] & groupbit) {

      qi = q[i];
      if (qi < qmin[type[i]]) {
        Hdia_inv[i] = 1. / (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1]);
        b_s[i]      = -((chi[type[i]]-2*qmin[type[i]]*omega[type[i]]) + chizj[i]);
      } else if (qi < qmax[type[i]]) {
        Hdia_inv[i] = 1. / (eta[type[i]]-s2d_self[type[i]-1]);
        b_s[i]      = -(chi[type[i]] + chizj[i]);
      } else {
        Hdia_inv[i] = 1. / (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1]);
        b_s[i]      = -((chi[type[i]]-2*qmax[type[i]]*omega[type[i]]) + chizj[i]);
      }

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
   threaded build of the ctip H matrix (full neighbor list).  Two-pass
   (count -> serial sum-scan offsets -> fill) reproduces the serial
   FixQEqCTIP::compute_H() packed CSR layout.
------------------------------------------------------------------------- */

void FixQEqCTIPOMP::compute_H_thr()
{
  int *ilist, *numneigh, **firstneigh;
  double **x = atom->x;
  int *mask = atom->mask;
  int *type = atom->type;

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
      const int itype = type[i];
      int mfill = H.firstnbr[i];
      for (int jj = 0; jj < jnum; jj++) {
        const int j = jlist[jj] & NEIGHMASK;
        const double dx = x[j][0] - x[i][0];
        const double dy = x[j][1] - x[i][1];
        const double dz = x[j][2] - x[i][2];
        const double r_sqr = dx*dx + dy*dy + dz*dz;
        if (r_sqr <= cutoff_sq) {
          const int jtype = type[j];
          const double r = sqrt(r_sqr);
          const double reff = cbrt(r_sqr * r + 1 / shieldcu[itype-1][jtype-1]);
          const double erfcd = exp(-cdamp * cdamp * r_sqr);
          const double tt = 1.0 / (1.0 + EWALD_P * cdamp * r);
          const double erfcc = tt * (A1 + tt * (A2 + tt * (A3 + tt * (A4 + tt * A5)))) * erfcd;
          H.jlist[mfill] = j;
          H.val[mfill] = 0.5*force->qqr2e*(erfcc/r+1/reff-1/r-e_shift[itype-1][jtype-1]
                          +f_shift[itype-1][jtype-1]*(r-cutoff)
                          -s2d_shift[itype-1][jtype-1]*0.5*(r-cutoff)*(r-cutoff));
          mfill++;
        }
      }
    }
  }

  if (m_fill >= H.m)
    error->all(FLERR, Error::NOLASTLINE,
               "Fix qeq/ctip/omp has insufficient H matrix size: m_fill={} H.m={}\n",m_fill, H.m);
}

/* ----------------------------------------------------------------------
   threaded sparse matrix-vector product (ctip qmin/qmax diagonal)
------------------------------------------------------------------------- */

void FixQEqCTIPOMP::sparse_matvec(sparse_matrix *A, double *x, double *b)
{
  const int nthreads = comm->nthreads;
  double *q = atom->q;
  int *type = atom->type;

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
    for (int i = 0; i < nloc; ++i) {
      if (atom->mask[i] & groupbit) {
        const double qi = q[i];
        if (qi < qmin[type[i]]) {
          b[i] = (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1])*x[i];
        } else if (qi < qmax[type[i]]) {
          b[i] = (eta[type[i]]-s2d_self[type[i]-1]) * x[i];
        } else {
          b[i] = (eta[type[i]]+2*omega[type[i]]-s2d_self[type[i]-1])*x[i];
        }
      }
    }

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
