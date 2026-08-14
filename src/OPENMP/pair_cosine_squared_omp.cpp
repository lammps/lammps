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

#include "pair_cosine_squared_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCosineSquaredOMP::PairCosineSquaredOMP(LAMMPS *lmp) :
    PairCosineSquared(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCosineSquaredOMP::compute(int eflag, int vflag)
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
void PairCosineSquaredOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *_noalias const special_lj = force->special_lj;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double fxtmp, fytmp, fztmp;

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const double r = sqrt(rsq);
        double evdwl = 0.0;
        double fpair;

        if (r <= sigma[itype][jtype]) {
          if (wcaflag[itype][jtype]) {
            const double r2inv = 1.0 / rsq;
            const double r6inv = r2inv * r2inv * r2inv;
            const double force_lj = r6inv * (lj12_f[itype][jtype] * r6inv - lj6_f[itype][jtype]);
            fpair = factor_lj * force_lj * r2inv;
            if (EFLAG) {
              evdwl = factor_lj * r6inv * (lj12_e[itype][jtype] * r6inv - lj6_e[itype][jtype]);
              if (sigma[itype][jtype] == cut[itype][jtype]) evdwl += factor_lj * epsilon[itype][jtype];
            }
          } else {
            fpair = 0.0;
            if (EFLAG) evdwl = -factor_lj * epsilon[itype][jtype];
          }
        } else {
          const double force_cos = -(MY_PI * epsilon[itype][jtype] / (2.0 * w[itype][jtype])) *
              sin(MY_PI * (r - sigma[itype][jtype]) / w[itype][jtype]);
          fpair = factor_lj * force_cos / r;
          if (EFLAG) {
            const double cosone = cos(MY_PI * (r - sigma[itype][jtype]) / (2.0 * w[itype][jtype]));
            evdwl = -factor_lj * epsilon[itype][jtype] * cosone * cosone;
          }
        }

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairCosineSquaredOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCosineSquared::memory_usage();
  return bytes;
}
