// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_morse_soft_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMorseSoftOMP::PairMorseSoftOMP(LAMMPS *lmp) :
  PairMorseSoft(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairMorseSoftOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (evflag) {
      if (eflag) {
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairMorseSoftOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  int i, j, ii, jj, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r, dr, dexp, dexp2, dexp3, factor_lj;
  double ea, phi, V0, iea2;
  double D, a, x0, l, B, s1, llf;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_lj = force->special_lj;

  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        dr = r - r0[itype][jtype];

        D = d0[itype][jtype];
        a = alpha[itype][jtype];
        x0 = r0[itype][jtype];
        dexp = exp(-a * dr);
        dexp2 = dexp * dexp;
        dexp3 = dexp2 * dexp;

        l = lambda[itype][jtype];

        ea = exp(a * x0);
        iea2 = exp(-2. * a * x0);

        V0 = D * dexp * (dexp - 2.0);
        B = -2.0 * D * iea2 * (ea - 1.0) / 3.0;

        if (l >= shift_range) {
          s1 = (l - 1.0) / (shift_range - 1.0);
          phi = V0 + B * dexp3 * s1;

          // Force computation:
          fpair = 3.0 * a * B * dexp3 * s1 + 2.0 * a * D * (dexp2 - dexp);
          fpair /= r;
        } else {
          llf = MathSpecial::powint(l / shift_range, nlambda);
          phi = V0 + B * dexp3;
          phi *= llf;

          // Force computation:
          if (r == 0.0) {
            fpair = 0.0;
          } else {
            fpair = 3.0 * a * B * dexp3 + 2.0 * a * D * (dexp2 - dexp);
            fpair *= llf / r;
          }
        }

        fpair *= factor_lj;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;

        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (EFLAG) { evdwl = phi * factor_lj; }

        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, 0.0, fpair, delx, dely, delz, thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairMorseSoftOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairMorseSoft::memory_usage();

  return bytes;
}
