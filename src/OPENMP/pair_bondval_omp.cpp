/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_bondval_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairBondValOMP::PairBondValOMP(LAMMPS *lmp) : PairBondVal(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  maxthreads = 0;
}

/* ---------------------------------------------------------------------- */

void PairBondValOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow the bond valence sum and its derivative if necessary
  // the sum holds one copy per thread, each atom->nmax long

  if ((atom->nmax > nmax) || (nthreads > maxthreads)) {
    memory->destroy(s0);
    memory->destroy(fp);
    nmax = atom->nmax;
    maxthreads = nthreads;
    memory->create(s0, (bigint) nthreads * nmax, "pair:s0");
    memory->create(fp, nmax, "pair:fp");
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

    if (evflag) {
      if (eflag) {
        if (force->newton_pair)
          eval<1, 1, 1>(ifrom, ito, thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair)
          eval<1, 0, 1>(ifrom, ito, thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair)
        eval<0, 0, 1>(ifrom, ito, thr);
      else
        eval<0, 0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairBondValOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int tid = thr->get_tid();

  // zero out this thread's copy of the bond valence sum

  const int nlen = NEWTON_PAIR ? nall : nlocal;
  double *_noalias const s0_t = s0 + (bigint) tid * nlen;
  memset(s0_t, 0, nlen * sizeof(double));

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      const int jtype = type[j];
      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;

      if (alpha[itype][jtype] != 0.0) {
        if (rsq < cutsq[itype][jtype]) {
          const double r = sqrt(rsq);
          s0_t[i] += pow(r0[itype][jtype] / r, alpha[itype][jtype]);
          if (NEWTON_PAIR || j < nlocal) s0_t[j] += pow(r0[jtype][itype] / r, alpha[jtype][itype]);
        }
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // reduce per thread bond valence sums, then communicate and sum them

  thr->timer(Timer::PAIR);
  data_reduce_thr(s0, nlen, nthreads, 1, tid);

  // wait until reduction is complete
  sync_threads();

  if (NEWTON_PAIR) {
#if defined(_OPENMP)
#if _OPENMP >= 202011
#pragma omp masked
#else
#pragma omp master
#endif
#endif
    { comm->reverse_comm(this); }

    // wait until master thread is done with communication
    sync_threads();
  }

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const double s = s0[i] - v0[itype][itype];
    fp[i] = sparam[itype][itype] * power_global * s;
    if (EFLAG) {
      const double phi = sparam[itype][itype] * s * s + energy0[itype];
      e_tally_thr(this, i, i, nlocal, NEWTON_PAIR, phi, 0.0, thr);
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // communicate derivative of the bond valence energy
  // MPI communication only on master thread
#if defined(_OPENMP)
#if _OPENMP >= 202011
#pragma omp masked
#else
#pragma omp master
#endif
#endif
  { comm->forward_comm(this); }

  // wait until master thread is done with communication
  sync_threads();

  // compute forces on each atom
  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp, fytmp, fztmp;
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;

      const int jtype = type[j];
      if (alpha[itype][jtype] != 0.0) {
        if (rsq < cutsq[itype][jtype]) {
          const double r = sqrt(rsq);
          const double recip = 1.0 / r;
          const double Aij =
              alpha[itype][jtype] * pow(r0[itype][jtype] * recip, alpha[itype][jtype]) * recip;
          const double psip = (fp[i] + fp[j]) * Aij;
          const double fpair = psip * recip;

          fxtmp += delx * fpair;
          fytmp += dely * fpair;
          fztmp += delz * fpair;
          if (NEWTON_PAIR || j < nlocal) {
            f[j].x -= delx * fpair;
            f[j].y -= dely * fpair;
            f[j].z -= delz * fpair;
          }

          // the pair terms carry no energy, it is all in the bond valence sums

          if (EVFLAG)
            ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, 0.0, fpair, delx, dely, delz, thr);
        }
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairBondValOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBondVal::memory_usage();
  bytes += (double) (maxthreads + 1) * nmax * sizeof(double);    // s0 copies and fp

  return bytes;
}
