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

#include "pair_coul_streitz_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulStreitzOMP::PairCoulStreitzOMP(LAMMPS *lmp) : PairCoulStreitz(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairCoulStreitzOMP::compute(int eflag, int vflag)
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
void PairCoulStreitzOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const int *_noalias const type = atom->type;
  const double *_noalias const q = atom->q;
  const double *_noalias const special_coul = force->special_coul;
  const int nlocal = atom->nlocal;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double fxtmp, fytmp, fztmp;

  if (kspacetype == 1) {

    // Wolf sum

    for (int ii = iifrom; ii < iito; ++ii) {
      const int i = ilist[ii];
      const double xtmp = x[i].x;
      const double ytmp = x[i].y;
      const double ztmp = x[i].z;
      const int itype = map[type[i]];
      const int iparam_i = elem1param[itype];
      const double qi = q[i];
      const double zei = params[iparam_i].zeta;
      fxtmp = fytmp = fztmp = 0.0;

      // self energy: ionization + wolf sum

      const double selfion = self(&params[iparam_i], qi);
      if (EFLAG)
        ev_tally_thr(this, i, i, nlocal, 0, 0.0, selfion, 0.0, 0.0, 0.0, 0.0, thr);

      // two-body interaction

      const int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;
        const int jtype = map[type[j]];
        const int iparam_j = elem1param[jtype];
        const double qj = q[j];
        const double zej = params[iparam_j].zeta;
        const double zj = params[iparam_j].zcore;

        double delr[3];
        delr[0] = xtmp - x[j].x;
        delr[1] = ytmp - x[j].y;
        delr[2] = ztmp - x[j].z;
        const double rsq = delr[0] * delr[0] + delr[1] * delr[1] + delr[2] * delr[2];

        if (rsq > cut_coulsq) continue;

        const double r = sqrt(rsq);

        double ci_jfi, dci_jfi, ci_fifj, dci_fifj;
        coulomb_integral_wolf(zei, zej, r, ci_jfi, dci_jfi, ci_fifj, dci_fifj);

        double ecoul, forcecoul;
        if (dsfflag == 1)
          fennell_sum(qi, qj, zj, r, ci_jfi, dci_jfi, ci_fifj, dci_fifj, ecoul, forcecoul);
        else
          wolf_sum(qi, qj, zj, r, ci_jfi, dci_jfi, ci_fifj, dci_fifj, ecoul, forcecoul);

        const double fpair = -forcecoul / r;

        fxtmp += delr[0] * fpair;
        fytmp += delr[1] * fpair;
        fztmp += delr[2] * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delr[0] * fpair;
          f[j].y -= delr[1] * fpair;
          f[j].z -= delr[2] * fpair;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, ecoul, fpair, delr[0], delr[1],
                       delr[2], thr);
      }
      f[i].x += fxtmp;
      f[i].y += fytmp;
      f[i].z += fztmp;
    }

  } else if (kspacetype == 2) {

    // Ewald sum

    for (int ii = iifrom; ii < iito; ++ii) {
      const int i = ilist[ii];
      const double xtmp = x[i].x;
      const double ytmp = x[i].y;
      const double ztmp = x[i].z;
      const int itype = map[type[i]];
      const int iparam_i = elem1param[itype];
      const double qi = q[i];
      const double zei = params[iparam_i].zeta;
      fxtmp = fytmp = fztmp = 0.0;

      // self ionization energy

      const double selfion = self(&params[iparam_i], qi);
      if (EFLAG)
        ev_tally_thr(this, i, i, nlocal, 0, 0.0, selfion, 0.0, 0.0, 0.0, 0.0, thr);

      // two-body interaction

      const int *jlist = firstneigh[i];
      const int jnum = numneigh[i];

      for (int jj = 0; jj < jnum; jj++) {
        int j = jlist[jj];
        j &= NEIGHMASK;
        const int jtype = map[type[j]];
        const int iparam_j = elem1param[jtype];
        const double qj = q[j];
        const double zej = params[iparam_j].zeta;
        const double zj = params[iparam_j].zcore;
        const double factor_coul = special_coul[sbmask(j)];

        double delr[3];
        delr[0] = xtmp - x[j].x;
        delr[1] = ytmp - x[j].y;
        delr[2] = ztmp - x[j].z;
        const double rsq = delr[0] * delr[0] + delr[1] * delr[1] + delr[2] * delr[2];

        if (rsq > cut_coulsq) continue;

        const double r = sqrt(rsq);

        double ci_jfi, dci_jfi, ci_fifj, dci_fifj;
        coulomb_integral_ewald(zei, zej, r, ci_jfi, dci_jfi, ci_fifj, dci_fifj);

        double ecoul, forcecoul;
        ewald_sum(qi, qj, zj, r, ci_jfi, dci_jfi, ci_fifj, dci_fifj, ecoul, forcecoul,
                  factor_coul);

        const double fpair = -forcecoul / r;

        fxtmp += delr[0] * fpair;
        fytmp += delr[1] * fpair;
        fztmp += delr[2] * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j].x -= delr[0] * fpair;
          f[j].y -= delr[1] * fpair;
          f[j].z -= delr[2] * fpair;
        }

        if (EVFLAG)
          ev_tally_thr(this, i, j, nlocal, NEWTON_PAIR, 0.0, ecoul, fpair, delr[0], delr[1],
                       delr[2], thr);
      }
      f[i].x += fxtmp;
      f[i].y += fytmp;
      f[i].z += fztmp;
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairCoulStreitzOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairCoulStreitz::memory_usage();
  return bytes;
}
