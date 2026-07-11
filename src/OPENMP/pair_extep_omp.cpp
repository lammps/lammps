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

#include "pair_extep_omp.h"

#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairExTePOMP::PairExTePOMP(LAMMPS *lmp) : PairExTeP(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairExTePOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  // SR_neigh must be called before the OMP parallel region
  SR_neigh();

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
void PairExTePOMP::eval(int iifrom, int iito, ThrData *const thr)
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

    // two-body interactions, skip half of them

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
      double delx = xtmp - x[j].x;
      double dely = ytmp - x[j].y;
      double delz = ztmp - x[j].z;
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

      double zeta_ij = 0.0;

      /* F_IJ (1) */
      double FXY, dFXY_dNdij, dFXY_dNdji, fa, fa_d, deng, fpair;
      double Ntij = Nt[i];
      double Ndij = Nd[i];
      double Ntji = Nt[j];
      double Ndji = Nd[j];
      const double r = sqrt(rsq1);
      const double fc_ij = ters_fc(r, &params[iparam_ij]);

      Ntij -= fc_ij;
      Ntji -= fc_ij;
      if (jtype != itype) {
        Ndij -= fc_ij;
        Ndji -= fc_ij;
      }
      if (Ntij < 0.0) Ntij = 0.0;
      if (Ndij < 0.0) Ndij = 0.0;
      if (Ntji < 0.0) Ntji = 0.0;
      if (Ndji < 0.0) Ndji = 0.0;
      FXY = F_corr(itype, jtype, Ndij, Ndji, &dFXY_dNdij, &dFXY_dNdji);

      double fenv, dfenv_ij;
      fenv = envelop_function(Ntij, Ntji, &dfenv_ij);

      const double Fc = fenv * FXY;
      const double dFc_dNtij = dfenv_ij * FXY;
      const double dFc_dNdij = fenv * dFXY_dNdij;

      fa = ters_fa(r, &params[iparam_ij]);
      fa_d = ters_fa_d(r, &params[iparam_ij]);
      deng = 0.5 * fa * Fc;
      fpair = 0.5 * fa_d * Fc / r;

      fxtmp += delr1[0] * fpair;
      fytmp += delr1[1] * fpair;
      fztmp += delr1[2] * fpair;
      f[j].x -= delr1[0] * fpair;
      f[j].y -= delr1[1] * fpair;
      f[j].z -= delr1[2] * fpair;

      if (EVFLAG)
        ev_tally_thr(this, i, j, nlocal, /* newton_pair */ 1, deng, 0.0, -fpair, -delr1[0],
                     -delr1[1], -delr1[2], thr);
      /* END F_IJ (1) */

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

        const double r2 = sqrt(rsq2);

        zeta_ij += zeta(&params[iparam_ijk], r, r2, delr1, delr2);

        /* F_IJ (2) */
        const int iparam_ik = elem3param[itype][ktype][0];
        const double fc_ik_d = ters_fc_d(r2, &params[iparam_ik]);
        const double fc_prefac_ik_0 = fc_ik_d * fa / r2;
        double fc_prefac_ik = dFc_dNtij * fc_prefac_ik_0;
        fxtmp += fc_prefac_ik * delr2[0];
        fytmp += fc_prefac_ik * delr2[1];
        fztmp += fc_prefac_ik * delr2[2];
        f[k].x -= fc_prefac_ik * delr2[0];
        f[k].y -= fc_prefac_ik * delr2[1];
        f[k].z -= fc_prefac_ik * delr2[2];
        if (VFLAG_EITHER) v_tally2_thr(this, i, k, -fc_prefac_ik, delr2, thr);
        if (itype != ktype) {
          fc_prefac_ik = dFc_dNdij * fc_prefac_ik_0;
          fxtmp += fc_prefac_ik * delr2[0];
          fytmp += fc_prefac_ik * delr2[1];
          fztmp += fc_prefac_ik * delr2[2];
          f[k].x -= fc_prefac_ik * delr2[0];
          f[k].y -= fc_prefac_ik * delr2[1];
          f[k].z -= fc_prefac_ik * delr2[2];
          if (VFLAG_EITHER) v_tally2_thr(this, i, k, -fc_prefac_ik, delr2, thr);
        }
        /* END F_IJ (2) */
      }

      // pairwise force due to zeta

      double zprefactor = 0.0;
      double evdwl = 0.0;
      force_zeta(&params[iparam_ij], r, zeta_ij, fpair, zprefactor, EFLAG, evdwl);

      fxtmp += delr1[0] * fpair;
      fytmp += delr1[1] * fpair;
      fztmp += delr1[2] * fpair;
      f[j].x -= delr1[0] * fpair;
      f[j].y -= delr1[1] * fpair;
      f[j].z -= delr1[2] * fpair;

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
        attractive(&params[iparam_ijk], zprefactor, rsq1, rsq2, delr1, delr2, fi, fj, fk);

        fxtmp += fi[0];
        fytmp += fi[1];
        fztmp += fi[2];
        f[j].x += fj[0];
        f[j].y += fj[1];
        f[j].z += fj[2];
        f[k].x += fk[0];
        f[k].y += fk[1];
        f[k].z += fk[2];

        if (VFLAG_EITHER) v_tally3_thr(this, i, j, k, fj, fk, delr1, delr2, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairExTePOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairExTeP::memory_usage();
  return bytes;
}
