/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_lj_cut_sphere_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include "omp_compat.h"
using namespace LAMMPS_NS;
using MathSpecial::powint;
using MathSpecial::square;

/* ---------------------------------------------------------------------- */

PairLJCutSphereOMP::PairLJCutSphereOMP(LAMMPS *lmp) : PairLJCutSphere(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void PairLJCutSphereOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

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

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCutSphereOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const double *_noalias const radius = atom->radius;
  const int *_noalias const type = atom->type;
  const double *_noalias const special_lj = force->special_lj;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double xtmp, ytmp, ztmp, rtmp, delx, dely, delz, fxtmp, fytmp, fztmp;
  double rcutsq, rsq, r2inv, r6inv, forcelj, factor_lj, evdwl, sigma, sigma6, fpair;

  const int nlocal = atom->nlocal;
  int j, jj, jnum, jtype;

  evdwl = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const double *_noalias const epsiloni = epsilon[itype];
    const double *_noalias const cutsqi = cutsq[itype];

    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    rtmp = radius[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {

        // cutsq is maximum cutoff per type. Now compute and apply real cutoff

        sigma = 2.0 * mix_distance(rtmp, radius[j]);
        rcutsq = square(cut[itype][jtype] * sigma);

        if (rsq < rcutsq) {

          r2inv = 1.0 / rsq;
          r6inv = r2inv * r2inv * r2inv;

          sigma6 = powint(sigma, 6);
          forcelj = r6inv * 24.0 * epsiloni[jtype] * (2.0 * sigma6 * sigma6 * r6inv - sigma6);
          fpair = factor_lj * forcelj * r2inv;

          fxtmp += delx * fpair;
          fytmp += dely * fpair;
          fztmp += delz * fpair;
          if (NEWTON_PAIR || j < nlocal) {
            f[j].x -= delx * fpair;
            f[j].y -= dely * fpair;
            f[j].z -= delz * fpair;
          }

          if (EFLAG) {
            evdwl = r6inv * 4.0 * epsiloni[jtype];
            evdwl *= sigma6 * sigma6 * r6inv - sigma6;
            if (offset_flag && (cutsqi[jtype] > 0.0)) {
              const double ratio6 = sigma6 / powint(rcutsq, 3);
              evdwl -= 4.0 * epsiloni[jtype] * (ratio6 * ratio6 - ratio6);
            }
            evdwl *= factor_lj;
          }

          if (EVFLAG)
            ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair, delx, dely, delz, thr);
        }
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutSphereOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJCutSphere::memory_usage();

  return bytes;
}
