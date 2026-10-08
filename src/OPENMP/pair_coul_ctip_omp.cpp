// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Chun-Yeol You (DGIST)
------------------------------------------------------------------------- */

#include "pair_coul_ctip_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;

static constexpr double MY_PIS = 1.7724538509055159;
static constexpr double EWALD_P = 0.3275911;
static constexpr double A1 =  0.254829592;
static constexpr double A2 = -0.284496736;
static constexpr double A3 =  1.421413741;
static constexpr double A4 = -1.453152027;
static constexpr double A5 =  1.061405429;

/* ---------------------------------------------------------------------- */

PairCoulCTIPOMP::PairCoulCTIPOMP(LAMMPS *lmp) :
  PairCoulCTIP(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCoulCTIPOMP::compute(int eflag, int vflag)
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
        if (force->newton_pair) eval<1,1,1>(ifrom, ito, thr);
        else                    eval<1,1,0>(ifrom, ito, thr);
      } else {
        if (force->newton_pair) eval<1,0,1>(ifrom, ito, thr);
        else                    eval<1,0,0>(ifrom, ito, thr);
      }
    } else {
      if (force->newton_pair) eval<0,0,1>(ifrom, ito, thr);
      else                    eval<0,0,0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairCoulCTIPOMP::eval(int iifrom, int iito, ThrData * const thr)
{
  const auto * _noalias const x = (dbl3_t *) atom->x[0];
  auto * _noalias const f = (dbl3_t *) thr->get_f()[0];
  const double * _noalias const q = atom->q;
  const int * _noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * _noalias const special_coul = force->special_coul;
  const double qqrd2e = force->qqrd2e;

  const double erfcd_cut = exp(-alpha * alpha * cut_coulsq);
  const double t_cut = 1.0 / (1.0 + EWALD_P * alpha * cut_coul);
  const double erfcc_cut = t_cut * (A1 + t_cut * (A2 + t_cut * (A3 + t_cut * (A4 + t_cut * A5)))) * erfcd_cut;

  const int * _noalias const ilist    = list->ilist;
  const int * _noalias const numneigh = list->numneigh;
  int ** const firstneigh             = list->firstneigh;

  for (int ii = iifrom; ii < iito; ii++) {
    const int i = ilist[ii];
    const double qtmp = q[i];
    const double xtmp  = x[i].x;
    const double ytmp  = x[i].y;
    const double ztmp  = x[i].z;
    const int itype    = map[type[i]];
    const int iparam_i = elem1param[itype];
    const int * _noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp = 0.0, fytmp = 0.0, fztmp = 0.0;

    const double * _noalias const shieldcu_i  = shieldcu[iparam_i];
    const double * _noalias const reffc_i     = reffc[iparam_i];
    const double * _noalias const f_shift_i   = f_shift[iparam_i];
    const double * _noalias const s2d_shift_i = s2d_shift[iparam_i];

    const double selfion = self(&params[iparam_i], qtmp);

    if (EVFLAG) ev_tally_thr(this, i, i, nlocal, 0, 0.0, selfion, 0.0, 0.0, 0.0, 0.0, thr);

    if (EFLAG) {
      const double e_self = self_factor[iparam_i][iparam_i] * qtmp * qtmp;
      ev_tally_thr(this, i, i, nlocal, 0, 0.0, e_self, 0.0, 0.0, 0.0, 0.0, thr);
    }

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      const double factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq  = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_coulsq) {
        const int jtype    = map[type[j]];
        const int jparam_j = elem1param[jtype];
        const double r     = sqrt(rsq);
        const double reff  = cbrt(rsq * r + 1 / shieldcu_i[jparam_j]);
        const double reffsq = reff * reff;
        const double reff4  = reffsq * reffsq;
        const double prefactor = qqrd2e * qtmp * q[j] / r;
        const double erfcd = exp(-alpha * alpha * rsq);
        const double t     = 1.0 / (1.0 + EWALD_P * alpha * r);
        const double erfcc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * erfcd;

        double forcecoul = prefactor *
            (erfcc / r + 2.0 * alpha / MY_PIS * erfcd + rsq * r / reff4 - 1 / r -
             r * f_shift_i[jparam_j] + r * s2d_shift_i[jparam_j] * (r - cut_coul)) * r;
        if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
        const double fpair = forcecoul / rsq;

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
          ecoul = prefactor *
              (erfcc + r / reff - 1 - r * erfcc_cut / cut_coul - r / reffc_i[jparam_j] +
               r / cut_coul + r * f_shift_i[jparam_j] * (r - cut_coul) -
               r * s2d_shift_i[jparam_j] * (r - cut_coul) * (r - cut_coul) * 0.5);
          if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
        }

        if (EVFLAG) ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR,
                                  0.0, ecoul, fpair, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairCoulCTIPOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCoulCTIP::memory_usage();
  return bytes;
}
