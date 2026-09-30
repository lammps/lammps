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

#include "pair_coul_shield_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "interlayer_taper.h"
#include "math_special.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace InterLayer;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

PairCoulShieldOMP::PairCoulShieldOMP(LAMMPS *lmp) : PairCoulShield(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCoulShieldOMP::compute(int eflag, int vflag)
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
void PairCoulShieldOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const double *_noalias const q = atom->q;
  const tagint *_noalias const mol = atom->molecule;
  const int nlocal = atom->nlocal;
  const double *_noalias const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double fxtmp, fytmp, fztmp;

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double qtmp = q[i];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int itype = type[i];
    const tagint imol = mol[i];
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      if (rsq < cutsq[itype][jtype] && imol != mol[j]) {
        const double r = sqrt(rsq);
        const double r3 = rsq * r;
        const double rarg = 1.0 / sigmae[itype][jtype];
        const double th = r3 + cube(rarg);
        const double epsr = 1.0 / pow(th, 1.0 / 3.0);
        double depsdr = square(epsr);
        depsdr *= depsdr;
        const double Vc = qqrd2e * qtmp * q[j] * epsr;

        double Tap, dTap;
        if (tap_flag) {
          Tap = calc_Tap(r, cut[itype][jtype]);
          dTap = calc_dTap(r, cut[itype][jtype]);
        } else {
          Tap = 1.0;
          dTap = 0.0;
        }

        const double forcecoul = qqrd2e * qtmp * q[j] * r * depsdr;
        const double fvc = forcecoul * Tap - Vc * dTap / r;
        const double fpair = factor_coul * fvc;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        double ecoul = 0.0;
        if (EFLAG) {
          if (tap_flag)
            ecoul = Vc * Tap;
          else
            ecoul = Vc - offset[itype][jtype];
          ecoul *= factor_coul;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, ecoul, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairCoulShieldOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCoulShield::memory_usage();
  return bytes;
}
