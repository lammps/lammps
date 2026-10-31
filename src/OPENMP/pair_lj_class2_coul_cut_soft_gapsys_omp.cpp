/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

#include "pair_lj_class2_coul_cut_soft_gapsys_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJClass2CoulCutSoftGapsysOMP::PairLJClass2CoulCutSoftGapsysOMP(LAMMPS *lmp) :
    PairLJClass2CoulCutSoftGapsys(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJClass2CoulCutSoftGapsysOMP::compute(int eflag, int vflag)
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

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJClass2CoulCutSoftGapsysOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const double *_noalias const q = atom->q;
  const int *_noalias const type = atom->type;
  const double *_noalias const special_coul = force->special_coul;
  const double *_noalias const special_lj = force->special_lj;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int nlocal = atom->nlocal;
  const double qqrd2e = force->qqrd2e;
  double evdwl = 0.0;
  double ecoul = 0.0;

  // loop over neighbors of my atoms

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const double qtmp = q[i];
    const int itype = type[i];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    double fxtmp, fytmp, fztmp;
    fxtmp = fytmp = fztmp = 0.0;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_lj = special_lj[sbmask(j)];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        const double r2inv = 1.0 / rsq;
        const double cut_inner_q = (1.0 + sigmaq * fabs(qtmp * q[j])) * cut_factor[itype][jtype];

        // forcecoul and forcelj do not include the special bond factors.
        // those are applied once, when fpair, ecoul, and evdwl are assembled

        double forcecoul = 0.0;
        double forcelj = 0.0;
        if (EFLAG) ecoul = evdwl = 0.0;

        if (rsq < cut_inner_q * cut_inner_q) {
          const double r = sqrt(rsq);
          forcecoul = qqrd2e * qtmp * q[j];
          forcecoul *= rsq * (-2.0 / cut_inner_q + 3.0 / r) / (cut_inner_q * cut_inner_q);
          if (EFLAG) {
            ecoul = factor_coul * qqrd2e * qtmp * q[j];
            ecoul *= ((rsq / cut_inner_q - 3.0 * r) / cut_inner_q + 3.0) / cut_inner_q;
          }
        } else if (rsq < cut_coulsq[itype][jtype]) {
          forcecoul = qqrd2e * qtmp * q[j] * sqrt(r2inv);
          if (EFLAG) ecoul = factor_coul * qqrd2e * qtmp * q[j] * sqrt(r2inv);
        }

        if (rsq < cut_lj_inner[itype][jtype] * cut_lj_inner[itype][jtype]) {
          forcelj = rsq * (coeff_b1[itype][jtype] + coeff_b2[itype][jtype] / sqrt(rsq));
          if (EFLAG) {
            evdwl = -0.5 * coeff_b1[itype][jtype] * rsq - coeff_b2[itype][jtype] * sqrt(rsq) +
                coeff_a3[itype][jtype] - offset[itype][jtype];
            evdwl *= factor_lj;
          }
        } else if (rsq < cut_ljsq[itype][jtype]) {
          const double rinv = sqrt(r2inv);
          const double r3inv = r2inv * rinv;
          const double r6inv = r3inv * r3inv;
          forcelj = r6inv * (lj1[itype][jtype] * r3inv - lj2[itype][jtype]);
          if (EFLAG) {
            evdwl = r6inv * (lj3[itype][jtype] * r3inv - lj4[itype][jtype]) - offset[itype][jtype];
            evdwl *= factor_lj;
          }
        }

        const double fpair = (factor_coul * forcecoul + factor_lj * forcelj) * r2inv;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delx * fpair;
          f[j].y -= dely * fpair;
          f[j].z -= delz * fpair;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, ecoul, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJClass2CoulCutSoftGapsysOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJClass2CoulCutSoftGapsys::memory_usage();

  return bytes;
}
