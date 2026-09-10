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

#include "pair_lj_pirani_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using MathSpecial::square;

/* ---------------------------------------------------------------------- */

PairLJPiraniOMP::PairLJPiraniOMP(LAMMPS *lmp) : PairLJPirani(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJPiraniOMP::compute(int eflag, int vflag)
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
void PairLJPiraniOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const double *_noalias const special_lj = force->special_lj;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double xtmp, ytmp, ztmp, delx, dely, delz, fxtmp, fytmp, fztmp;
  double rsq, forceilj, factor_lj, evdwl, fpair;

  const int nlocal = atom->nlocal;
  int j, jj, jnum, jtype;

  evdwl = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const double *_noalias const cutsqi = cutsq[itype];
    const double *_noalias const offseti = offset[itype];
    const double *_noalias const alphai = alpha[itype];
    const double *_noalias const betai = beta[itype];
    const double *_noalias const gammai = gamma[itype];
    const double *_noalias const rmi = rm[itype];
    const double *_noalias const epsiloni = epsilon[itype];

    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
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
        const double r = sqrt(rsq);

        const double rx = r / rmi[jtype];
        const double n_x = alphai[jtype] * rx * rx + betai[jtype];
        const double pow_rx_n_x = pow(1.0 / rx, n_x);
        const double pow_rx_gamma = pow(1.0 / rx, gammai[jtype]);

        double filj1 = -2.0 * alphai[jtype] * gammai[jtype] * rx * pow_rx_n_x /
            (square(n_x - gammai[jtype]) * rmi[jtype]);

        double filj2 = 2.0 * alphai[jtype] * rx * n_x * pow_rx_gamma /
            (square(n_x - gammai[jtype]) * rmi[jtype]);

        double filj3 =
            -2.0 * alphai[jtype] * rx * pow_rx_gamma / (rmi[jtype] * (n_x - gammai[jtype]));

        double filj4 = 2.0 * alphai[jtype] * gammai[jtype] * (rx / rmi[jtype]) * log(1 / rx) *
            pow_rx_n_x / (n_x - gammai[jtype]);

        double filj5 = -1.0 * gammai[jtype] * n_x * pow_rx_n_x / (r * (n_x - gammai[jtype]));

        double filj6 = 1.0 * gammai[jtype] * n_x * pow_rx_gamma / (r * (n_x - gammai[jtype]));

        // F = -dV/dr
        forceilj = -epsiloni[jtype] * (filj1 + filj2 + filj3 + filj4 + filj5 + filj6);
        fpair = factor_lj * forceilj / r;    // F_x = -x/r * dV/dr (chain rule)

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        if (EFLAG) {
          double ilj1 = epsiloni[jtype] * gammai[jtype] * pow(1 / rx, n_x) / (n_x - gammai[jtype]);
          double ilj2 = -epsiloni[jtype] * n_x * pow(1 / rx, gammai[jtype]) / (n_x - gammai[jtype]);

          evdwl = ilj1 + ilj2 - offseti[jtype];
          evdwl *= factor_lj;
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

double PairLJPiraniOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJPirani::memory_usage();

  return bytes;
}
