/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_bondval_vec_omp.h"

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

PairBondValVecOMP::PairBondValVecOMP(LAMMPS *lmp) : PairBondValVec(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  maxthreads = 0;
}

/* ---------------------------------------------------------------------- */

void PairBondValVecOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

  // grow the bond valence vector sums and their derivatives if necessary
  // the sums hold one copy per thread, each atom->nmax long

  if ((atom->nmax > nmax) || (nthreads > maxthreads)) {
    memory->destroy(s0x);
    memory->destroy(s0y);
    memory->destroy(s0z);
    memory->destroy(Dix);
    memory->destroy(Diy);
    memory->destroy(Diz);
    nmax = atom->nmax;
    maxthreads = nthreads;
    memory->create(s0x, (bigint) nthreads * nmax, "pair:s0x");
    memory->create(s0y, (bigint) nthreads * nmax, "pair:s0y");
    memory->create(s0z, (bigint) nthreads * nmax, "pair:s0z");
    memory->create(Dix, nmax, "pair:Dix");
    memory->create(Diy, nmax, "pair:Diy");
    memory->create(Diz, nmax, "pair:Diz");
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
void PairBondValVecOMP::eval(int iifrom, int iito, ThrData *const thr)
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

  // zero out this thread's copy of the bond valence vector sum

  const int nlen = NEWTON_PAIR ? nall : nlocal;
  double *_noalias const s0x_t = s0x + (bigint) tid * nlen;
  double *_noalias const s0y_t = s0y + (bigint) tid * nlen;
  double *_noalias const s0z_t = s0z + (bigint) tid * nlen;
  memset(s0x_t, 0, nlen * sizeof(double));
  memset(s0y_t, 0, nlen * sizeof(double));
  memset(s0z_t, 0, nlen * sizeof(double));

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

      if (rsq < cutsq[itype][jtype]) {
        const double recip = 1.0 / sqrt(rsq);
        s0x_t[i] += pow(r0[itype][jtype] * recip, alpha[itype][jtype]) * recip * (delx);
        s0y_t[i] += pow(r0[itype][jtype] * recip, alpha[itype][jtype]) * recip * (dely);
        s0z_t[i] += pow(r0[itype][jtype] * recip, alpha[itype][jtype]) * recip * (delz);

        if (NEWTON_PAIR || j < nlocal) {
          s0x_t[j] -= pow(r0[jtype][itype] * recip, alpha[jtype][itype]) * recip * (delx);
          s0y_t[j] -= pow(r0[jtype][itype] * recip, alpha[jtype][itype]) * recip * (dely);
          s0z_t[j] -= pow(r0[jtype][itype] * recip, alpha[jtype][itype]) * recip * (delz);
        }
      }
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // reduce per thread bond valence vector sums, then communicate and sum them

  thr->timer(Timer::PAIR);
  data_reduce_thr(s0x, nlen, nthreads, 1, tid);
  data_reduce_thr(s0y, nlen, nthreads, 1, tid);
  data_reduce_thr(s0z, nlen, nthreads, 1, tid);

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
    const double s = s0x[i] * s0x[i] + s0y[i] * s0y[i] + s0z[i] * s0z[i] -
        bvvv0[itype][itype] * bvvv0[itype][itype];
    Dix[i] = bvvsparam[itype][itype] * power_global * 2.0 * s0x[i] * s;
    Diy[i] = bvvsparam[itype][itype] * power_global * 2.0 * s0y[i] * s;
    Diz[i] = bvvsparam[itype][itype] * power_global * 2.0 * s0z[i] * s;
    if (EFLAG) {
      const double phi = bvvsparam[itype][itype] * s * s;
      e_tally_thr(this, i, i, nlocal, NEWTON_PAIR, phi, 0.0, thr);
    }
  }

  // wait until all threads are done with computation
  sync_threads();

  // communicate derivatives of the bond valence vector energy
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
      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        const double recip = 1.0 / r;
        const double recip2 = recip * recip;
        const double Aij = pow(r0[itype][jtype] * recip, alpha[itype][jtype]) * recip;
        const double Eij = (alpha[itype][jtype] + 1.0) * recip2;

        const double fx = (Dix[j] - Dix[i]) * Aij + (Dix[i] - Dix[j]) * Eij * delx * delx * Aij +
            (Diy[i] - Diy[j]) * Eij * delx * dely * Aij +
            (Diz[i] - Diz[j]) * Eij * delx * delz * Aij;
        const double fy = (Diy[j] - Diy[i]) * Aij + (Diy[i] - Diy[j]) * Eij * dely * dely * Aij +
            (Diz[i] - Diz[j]) * Eij * dely * delz * Aij +
            (Dix[i] - Dix[j]) * Eij * dely * delx * Aij;
        const double fz = (Diz[j] - Diz[i]) * Aij + (Diz[i] - Diz[j]) * Eij * delz * delz * Aij +
            (Dix[i] - Dix[j]) * Eij * delz * delx * Aij +
            (Diy[i] - Diy[j]) * Eij * delz * dely * Aij;

        fxtmp += fx;
        fytmp += fy;
        fztmp += fz;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= fx;
          f[j].y -= fy;
          f[j].z -= fz;
        }

        // the pair terms carry no energy, it is all in the bond valence vector sums

        if (EVFLAG)
          ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, 0.0, fx, fy, fz, delx, dely, delz,
                           thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairBondValVecOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairBondValVec::memory_usage();
  bytes += (double) 3 * (maxthreads + 1) * nmax * sizeof(double);    // s0[xyz] copies, Di[xyz]

  return bytes;
}
