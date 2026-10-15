/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Claude Fable 5.1 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_lj_long_dipole_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace EwaldConst;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJLongDipoleLongOMP::PairLJLongDipoleLongOMP(LAMMPS *lmp) :
    PairLJLongDipoleLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJLongDipoleLongOMP::compute(int eflag, int vflag)
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
void PairLJLongDipoleLongOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const double *const *const x = atom->x;
  double *const *const f = thr->get_f();
  const double *const *const mu = atom->mu;
  double *const *const tq = thr->get_torque();
  const double *const q = atom->q;
  const int *const type = atom->type;
  const double *const special_coul = force->special_coul;
  const double *const special_lj = force->special_lj;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const int nlocal = atom->nlocal;
  const double qqrd2e = force->qqrd2e;
  const int order3 = ewald_order & EWALD_DIPOLE, order6 = ewald_order & EWALD_DISP;
  const double g2 = g_ewald * g_ewald, g6 = g2 * g2 * g2, g8 = g6 * g2;
  const double C1 = 2.0 * g_ewald / MY_PIS;
  const double C2 = 2.0 * g2 * C1;
  const double C3 = 2.0 * g2 * C2;

  double evdwl = 0.0, ecoul = 0.0;
  double force_coul, force_lj;
  double B0, B1, B2, B3, G0, G1, G2, mudi, mudj, muij;
  double mui[3], muj[3] = {0.0, 0.0, 0.0}, xi[3], d[3];

  for (int ii = iifrom; ii < iito; ++ii) {    // loop over all neighs
    const int i = ilist[ii];
    double *const tqi = tq[i];
    const double qi = q[i];    // initialize constants
    const int typei = type[i];
    const double *const offseti = offset[typei];
    const double *const lj1i = lj1[typei];
    const double *const lj2i = lj2[typei];
    const double *const lj3i = lj3[typei];
    const double *const lj4i = lj4[typei];
    const double *const cutsqi = cutsq[typei];
    const double *const cut_ljsqi = cut_ljsq[typei];
    xi[0] = x[i][0];
    xi[1] = x[i][1];
    xi[2] = x[i][2];
    mui[0] = mu[i][0];
    mui[1] = mu[i][1];
    mui[2] = mu[i][2];

    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {    // loop over neighbors
      int j = jlist[jj];
      const int ni = sbmask(j);    // special index
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];    // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      const double rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
      const int typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
      const double r2inv = 1.0 / rsq;

      double force_d[3] = {0.0, 0.0, 0.0}, ti[3] = {0.0, 0.0, 0.0}, tj[3] = {0.0, 0.0, 0.0};

      if (order3 && (rsq < cut_coulsq)) {    // dipole
        muj[0] = mu[j][0];
        muj[1] = mu[j][1];
        muj[2] = mu[j][2];
        {    // series real space
          const double r = sqrt(rsq);
          const double grij = g_ewald * r;
          double fexp = exp(-grij * grij) * qqrd2e;

          B0 = 1.0 / (1.0 + EWALD_P * grij);    // eqn 2.8
          B0 *= ((((A5 * B0 + A4) * B0 + A3) * B0 + A2) * B0 + A1) * fexp / r;
          B1 = (B0 + C1 * fexp) * r2inv;
          B2 = (3.0 * B1 + C2 * fexp) * r2inv;
          B3 = (5.0 * B2 + C3 * fexp) * r2inv;

          mudi = mui[0] * d[0] + mui[1] * d[1] + mui[2] * d[2];
          mudj = muj[0] * d[0] + muj[1] * d[1] + muj[2] * d[2];
          muij = mui[0] * muj[0] + mui[1] * muj[1] + mui[2] * muj[2];
          const double qj = q[j];
          G0 = qi * qj;    // eqn 2.10
          G1 = qi * mudj - qj * mudi + muij;
          G2 = -mudi * mudj;
          force_coul = G0 * B1 + G1 * B2 + G2 * B3;

          mudi *= B2;
          mudj *= B2;    // torque contribs
          ti[0] = mudj * d[0] + (qj * d[0] - muj[0]) * B1;
          ti[1] = mudj * d[1] + (qj * d[1] - muj[1]) * B1;
          ti[2] = mudj * d[2] + (qj * d[2] - muj[2]) * B1;

          if (NEWTON_PAIR || j < nlocal) {
            tj[0] = mudi * d[0] - (qi * d[0] + mui[0]) * B1;
            tj[1] = mudi * d[1] - (qi * d[1] + mui[1]) * B1;
            tj[2] = mudi * d[2] - (qi * d[2] + mui[2]) * B1;
          }

          if (EFLAG) ecoul = G0 * B0 + G1 * B1 + G2 * B2;
          if (ni > 0) {    // adj part, eqn 2.13
            fexp = qqrd2e * (1.0 - special_coul[ni]) / r;
            force_coul -= fexp * ((3.0 * G1 + 15.0 * G2 * r2inv) * r2inv + G0) * r2inv;
            if (EFLAG) ecoul -= fexp * ((G1 + 3.0 * G2 * r2inv) * r2inv + G0);
            B1 -= fexp * r2inv;
          }
          B0 = mudj + qj * B1;
          B3 = -qi * B1 + mudi;    // position independent
          if (ni > 0) B0 -= fexp * 3.0 * mudj * r2inv * r2inv / B2;
          if (ni > 0) B3 -= fexp * 3.0 * mudi * r2inv * r2inv / B2;
          force_d[0] = B0 * mui[0] + B3 * muj[0];    // force contribs
          force_d[1] = B0 * mui[1] + B3 * muj[1];
          force_d[2] = B0 * mui[2] + B3 * muj[2];
          if (ni > 0) {
            ti[0] -= fexp * (3.0 * mudj * r2inv * r2inv * d[0] / B2 + (qj * r2inv * d[0] - muj[0] * r2inv));
            ti[1] -= fexp * (3.0 * mudj * r2inv * r2inv * d[1] / B2 + (qj * r2inv * d[1] - muj[1] * r2inv));
            ti[2] -= fexp * (3.0 * mudj * r2inv * r2inv * d[2] / B2 + (qj * r2inv * d[2] - muj[2] * r2inv));
            if (NEWTON_PAIR || j < nlocal) {
              tj[0] -= fexp * (3.0 * mudi * r2inv * r2inv * d[0] / B2 - (qi * r2inv * d[0] + mui[0] * r2inv));
              tj[1] -= fexp * (3.0 * mudi * r2inv * r2inv * d[1] / B2 - (qi * r2inv * d[1] + mui[1] * r2inv));
              tj[2] -= fexp * (3.0 * mudi * r2inv * r2inv * d[2] / B2 - (qi * r2inv * d[2] + mui[2] * r2inv));
            }
          }
        }    // table real space
      } else {
        force_coul = ecoul = 0.0;
      }

      if (rsq < cut_ljsqi[typej]) {    // lj
        if (order6) {                  // long-range lj
          const double r6inv = r2inv * r2inv * r2inv;
          const double r12inv = r6inv * r6inv;
          const double gr2 = g2 * rsq, a2 = 1.0 / gr2;
          const double expterm = a2 * exp(-gr2) * lj4i[typej];    // damped 1/r^6 reciprocal term
          const double g6term = g6 * ((a2 + 1.0) * a2 + 0.5) * expterm;
          const double g8term = g8 * (((6.0 * a2 + 6.0) * a2 + 3.0) * a2 + 1.0) * expterm * rsq;
          if (ni == 0) {
            force_lj = r12inv * lj1i[typej] - g8term;
            if (EFLAG) evdwl = r12inv * lj3i[typej] - g6term;
          } else {    // special case
            const double factor = special_lj[ni], t = r6inv * (1.0 - factor);
            force_lj = factor * r12inv * lj1i[typej] - g8term + t * lj2i[typej];
            if (EFLAG) evdwl = factor * r12inv * lj3i[typej] - g6term + t * lj4i[typej];
          }
        } else {    // cut lj
          const double r6inv = r2inv * r2inv * r2inv;
          if (ni == 0) {
            force_lj = r6inv * (r6inv * lj1i[typej] - lj2i[typej]);
            if (EFLAG) evdwl = r6inv * (r6inv * lj3i[typej] - lj4i[typej]) - offseti[typej];
          } else {    // special case
            const double factor = special_lj[ni];
            force_lj = factor * r6inv * (r6inv * lj1i[typej] - lj2i[typej]);
            if (EFLAG)
              evdwl = factor * (r6inv * (r6inv * lj3i[typej] - lj4i[typej]) - offseti[typej]);
          }
        }
        force_lj *= r2inv;
      } else
        force_lj = evdwl = 0.0;

      const double fpair = force_coul + force_lj;    // force
      const double fx = d[0] * fpair + force_d[0];
      const double fy = d[1] * fpair + force_d[1];
      const double fz = d[2] * fpair + force_d[2];
      f[i][0] += fx;
      f[i][1] += fy;
      f[i][2] += fz;
      tqi[0] += mui[1] * ti[2] - mui[2] * ti[1];    // torque
      tqi[1] += mui[2] * ti[0] - mui[0] * ti[2];
      tqi[2] += mui[0] * ti[1] - mui[1] * ti[0];
      if (NEWTON_PAIR || j < nlocal) {
        f[j][0] -= fx;
        f[j][1] -= fy;
        f[j][2] -= fz;
        double *const tqj = tq[j];
        tqj[0] += muj[1] * tj[2] - muj[2] * tj[1];
        tqj[1] += muj[2] * tj[0] - muj[0] * tj[2];
        tqj[2] += muj[0] * tj[1] - muj[1] * tj[0];
      }

      if (EVFLAG)
        ev_tally_xyz_thr(this, i, j, nlocal, NEWTON_PAIR, evdwl, ecoul, fx, fy, fz, d[0], d[1],
                         d[2], thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLJLongDipoleLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJLongDipoleLong::memory_usage();

  return bytes;
}
