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

#include "pair_gw_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGWOMP::PairGWOMP(LAMMPS *lmp) : PairGW(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairGWOMP::compute(int eflag, int vflag)
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
        if (vflag_either)
          eval<1, 1, 1>(ifrom, ito, thr);
        else
          eval<1, 1, 0>(ifrom, ito, thr);
      } else {
        if (vflag_either)
          eval<1, 0, 1>(ifrom, ito, thr);
        else
          eval<1, 0, 0>(ifrom, ito, thr);
      }
    } else
      eval<0, 0, 0>(ifrom, ito, thr);

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

template <int EVFLAG, int EFLAG, int VFLAG_EITHER>
void PairGWOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const tagint *_noalias const tag = atom->tag;
  const int *_noalias const type = atom->type;
  const int nlocal = atom->nlocal;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  double fxtmp, fytmp, fztmp;

  for (int ii = iifrom; ii < iito; ++ii) {
    const int i = ilist[ii];
    const tagint itag = tag[i];
    const int itype = map[type[i]];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const int *jlist = firstneigh[i];
    const int jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them (using FULL list + tag filter)

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const tagint jtag = tag[j];

      if (itag > jtag) {
        if ((itag + jtag) % 2 == 0) continue;
      } else if (itag < jtag) {
        if ((itag + jtag) % 2 == 1) continue;
      } else {
        if (x[j].z < ztmp) continue;
        if (x[j].z == ztmp && x[j].y < ytmp) continue;
        if (x[j].z == ztmp && x[j].y == ytmp && x[j].x < xtmp) continue;
      }

      const int jtype = map[type[j]];
      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int iparam_ij = elem3param[itype][jtype][jtype];
      if (rsq > params[iparam_ij].cutsq) continue;

      double evdwl = 0.0, fpair = 0.0;
      repulsive(&params[iparam_ij], rsq, fpair, EFLAG, evdwl);

      fxtmp += delx * fpair;
      fytmp += dely * fpair;
      fztmp += delz * fpair;
      f[j].x -= delx * fpair;
      f[j].y -= dely * fpair;
      f[j].z -= delz * fpair;

      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1, evdwl, 0.0, fpair, delx, dely, delz,
                     thr);
    }

    // three-body interactions

    double fjxtmp, fjytmp, fjztmp;

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      const int jtype = map[type[j]];
      const int iparam_ij = elem3param[itype][jtype][jtype];

      double delr1[3];
      delr1[0] = x[j].x - xtmp;
      delr1[1] = x[j].y - ytmp;
      delr1[2] = x[j].z - ztmp;
      const double rsq1 = delr1[0] * delr1[0] + delr1[1] * delr1[1] + delr1[2] * delr1[2];
      if (rsq1 > params[iparam_ij].cutsq) continue;

      // accumulate bondorder zeta for each i-j interaction via loop over k
      // GW starts at 1.0, not 0.0

      double zeta_ij = 1.0;

      for (int kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        int k = jlist[kk];
        k &= NEIGHMASK;
        const int ktype = map[type[k]];
        const int iparam_ijk = elem3param[itype][jtype][ktype];

        double delr2[3];
        delr2[0] = x[k].x - xtmp;
        delr2[1] = x[k].y - ytmp;
        delr2[2] = x[k].z - ztmp;
        const double rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        zeta_ij += zeta(&params[iparam_ijk], rsq1, rsq2, delr1, delr2);
      }

      // pairwise force due to zeta

      double fpair = 0.0, prefactor = 0.0, evdwl = 0.0;
      force_zeta(&params[iparam_ij], rsq1, zeta_ij, fpair, prefactor, EFLAG, evdwl);

      fjxtmp = fjytmp = fjztmp = 0.0;
      fxtmp += delr1[0] * fpair;
      fytmp += delr1[1] * fpair;
      fztmp += delr1[2] * fpair;
      fjxtmp -= delr1[0] * fpair;
      fjytmp -= delr1[1] * fpair;
      fjztmp -= delr1[2] * fpair;

      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1, evdwl, 0.0, -fpair, -delr1[0],
                     -delr1[1], -delr1[2], thr);

      // attractive term via loop over k

      for (int kk = 0; kk < jnum; kk++) {
        if (jj == kk) continue;
        int k = jlist[kk];
        k &= NEIGHMASK;
        const int ktype = map[type[k]];
        const int iparam_ijk = elem3param[itype][jtype][ktype];

        double delr2[3];
        delr2[0] = x[k].x - xtmp;
        delr2[1] = x[k].y - ytmp;
        delr2[2] = x[k].z - ztmp;
        const double rsq2 = delr2[0] * delr2[0] + delr2[1] * delr2[1] + delr2[2] * delr2[2];
        if (rsq2 > params[iparam_ijk].cutsq) continue;

        double fi[3], fj[3], fk[3];
        attractive(&params[iparam_ijk], prefactor, rsq1, rsq2, delr1, delr2, fi, fj, fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k].x += fk[0];
        f[k].y += fk[1];
        f[k].z += fk[2];

        if (VFLAG_EITHER) v_tally3_thr(this, i, j, k, fj, fk, delr1, delr2, thr);
      }
      f[j].x += fjxtmp;
      f[j].y += fjytmp;
      f[j].z += fjztmp;
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairGWOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairGW::memory_usage();
  return bytes;
}
