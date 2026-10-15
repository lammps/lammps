/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Claude Fable 5.1 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_smatb_single_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairSMATBSingleOMP::PairSMATBSingleOMP(LAMMPS *lmp) : PairSMATBSingle(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairSMATBSingleOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow the per-thread copies of the bonding energy array if necessary

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->grow(on_eb, nmax * nthreads, "pair_smatb:on_eb");
  }

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag, vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    // hook up and zero the per-thread copy of the bonding energy array
    thr->init_eam(nall, on_eb);

    if (evflag) {
      if (eflag) {
        eval<1, 1>(ifrom, ito, thr);
      } else {
        eval<1, 0>(ifrom, ito, thr);
      }
    } else {
      eval<0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG>
void PairSMATBSingleOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  double *const on_eb_t = thr->get_rho();
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int tid = thr->get_tid();
  const int nthreads = comm->nthreads;
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int newton_pair = force->newton_pair;

  double del[3], dij, polyval, polyval2, polyval3, polyval4, polyval5, qsiexpq;

  // FIRST LOOP: CALCULATES the squared bonding energy and accumulate it in on_eb for each atom
  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; ++jj) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      del[0] = xtmp - x[j].x;
      del[1] = ytmp - x[j].y;
      del[2] = ztmp - x[j].z;
      const double dijsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];

      if (dijsq < cutOffEnd2) {
        dij = sqrt(dijsq);
        if (dij < cutOffStart) {
          qsiexpq = (QSI * QSI) * exp(2.0 * q * (1.0 - dij / r0));
        } else {
          polyval = dij - cutOffEnd;
          polyval3 = polyval * polyval * polyval;
          polyval4 = polyval3 * polyval;
          polyval5 = polyval4 * polyval;
          qsiexpq = x5 * polyval5 + x4 * polyval4 + x3 * polyval3;
          qsiexpq = qsiexpq * qsiexpq;
        }
        on_eb_t[i] += qsiexpq;
        if (newton_pair) on_eb_t[j] += qsiexpq;
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // reduce the per-thread copies of the squared bonding energy
  thr->timer(Timer::PAIR);
  data_reduce_thr(on_eb, nall, nthreads, 1, tid);

  // wait until reduction is complete
  sync_threads();

  // communicate the squared bonding energy between the various bins
  // MPI communication only on master thread
  if (newton_pair) {
#if defined(_OPENMP)
#if _OPENMP >= 202011
#pragma omp masked
#else
#pragma omp master
#endif
#endif
    {
      comm->reverse_comm(this);
    }

    // wait until master thread is done with communication
    sync_threads();
  }

  // Support Loop: take the square root of the bonding energy and
  // accumulate it in the energy accumulator if needed the store the
  // reciprocal in on_eb in order to not do it in the SECOND LOOP

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    if (i < nlocal) {
      const double eb_i = sqrt(on_eb[i]);
      if (eb_i != 0.0) {
        on_eb[i] = 1.0 / eb_i;
      } else {
        on_eb[i] = 0.0;
      }
      // if needed the bonding energy is accumulated:
      if (EFLAG) e_tally_thr(this, i, i, nlocal, newton_pair, -eb_i, 0.0, thr);
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // this communication stores the denominators in the ghosts atoms,
  // this is needed because of how forces are calculated
  // MPI communication only on master thread
#if defined(_OPENMP)
#if _OPENMP >= 202011
#pragma omp masked
#else
#pragma omp master
#endif
#endif
  {
    comm->forward_comm(this);
  }

  // wait until master thread is done with communication
  sync_threads();

  // SECOND LOOP: given on_eb[i] calculates forces and energies
  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      del[0] = xtmp - x[j].x;
      del[1] = ytmp - x[j].y;
      del[2] = ztmp - x[j].z;

      const double dijsq = del[0] * del[0] + del[1] * del[1] + del[2] * del[2];
      if (dijsq < cutOffEnd2) {
        double aexpp, Fr, Fb;
        dij = sqrt(dijsq);
        if (dij < cutOffStart) {
          const double espo = 1.0 - dij / r0;
          aexpp = exp(p * espo) * A;
          Fr = (2.0 * aexpp) * (p / r0);
          qsiexpq = (QSI * QSI) * exp(2.0 * q * espo);
          Fb = -qsiexpq * q / r0;
        } else {
          polyval = dij - cutOffEnd;
          polyval2 = polyval * polyval;
          polyval3 = polyval2 * polyval;
          polyval4 = polyval3 * polyval;
          polyval5 = polyval4 * polyval;
          aexpp = a5 * polyval5 + a4 * polyval4 + a3 * polyval3;
          Fr = -2.0 * (5.0 * a5 * polyval4 + 4.0 * a4 * polyval3 + 3.0 * a3 * polyval2);
          qsiexpq = x5 * polyval5 + x4 * polyval4 + x3 * polyval3;
          Fb = (5.0 * x5 * polyval4 + 4.0 * x4 * polyval3 + 3.0 * x3 * polyval2) * qsiexpq;
        }

        // calculates the module of the pair energy between i and j
        const double fpair = (Fb * (on_eb[i] + on_eb[j]) + Fr) / dij;

        f[i].x += del[0] * fpair;
        f[i].y += del[1] * fpair;
        f[i].z += del[2] * fpair;
        if (newton_pair || j < nlocal) {
          f[j].x -= del[0] * fpair;
          f[j].y -= del[1] * fpair;
          f[j].z -= del[2] * fpair;
        }
        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, newton_pair, 2.0 * aexpp, 0.0, fpair, del[0], del[1],
                       del[2], thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairSMATBSingleOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairSMATBSingle::memory_usage();

  return bytes;
}
