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

#include "pair_lj_long_coul_long_dielectric_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace EwaldConst;
using MathConst::MY_PIS;

/* ---------------------------------------------------------------------- */

PairLJLongCoulLongDielectricOMP::PairLJLongCoulLongDielectricOMP(LAMMPS *_lmp) :
    PairLJLongCoulLongDielectric(_lmp), ThrOMP(_lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 0;
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLongDielectricOMP::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    memory->destroy(epot);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
    memory->create(epot, nmax, "pair:epot");
  }

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
        eval<1, 1>(ifrom, ito, thr);
      } else {
        eval<1, 0>(ifrom, ito, thr);
      }
    } else {
      eval<0, 0>(ifrom, ito, thr);
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  }    // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG>
void PairLJLongCoulLongDielectricOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  const auto *_noalias const x = (dbl3_t *) atom->x[0];
  auto *_noalias const f = (dbl3_t *) thr->get_f()[0];
  const double *_noalias const q = atom->q_scaled;
  const double *_noalias const eps = atom->epsilon;
  const auto *_noalias const norm = (dbl4_t *) atom->mu[0];
  const double *_noalias const curvature = atom->curvature;
  const double *_noalias const area = atom->area;
  const int *_noalias const type = atom->type;
  const double *_noalias const special_coul = force->special_coul;
  const double *_noalias const special_lj = force->special_lj;
  const int *_noalias const ilist = list->ilist;
  const int *_noalias const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  const double qqrd2e = force->qqrd2e;
  const int order1 = ewald_order & (1 << 1), order6 = ewald_order & (1 << 6);
  const double g2 = g_ewald_6 * g_ewald_6, g6 = g2 * g2 * g2, g8 = g6 * g2;

  int itable = 0;
  double evdwl = 0.0, ecoul = 0.0;
  double fraction = 0.0, table = 0.0;
  double grij, expm2, prefactor = 0.0, t, erfc = 0.0, prefactorE, efield_i, epot_i;
  double force_coul, force_lj;

  for (int ii = iifrom; ii < iito; ++ii) {    // loop over my atoms
    const int i = ilist[ii];
    const double qtmp = q[i];
    const double xtmp = x[i].x;
    const double ytmp = x[i].y;
    const double ztmp = x[i].z;
    const double etmp = eps[i];
    const int itype = type[i];

    const double *_noalias const offseti = offset[itype];
    const double *_noalias const lj1i = lj1[itype];
    const double *_noalias const lj2i = lj2[itype];
    const double *_noalias const lj3i = lj3[itype];
    const double *_noalias const lj4i = lj4[itype];
    const double *_noalias const cut_ljsqi = cut_ljsq[itype];
    const int *_noalias const jlist = firstneigh[i];
    const int jnum = numneigh[i];
    double fxtmp, fytmp, fztmp, extmp, eytmp, eztmp;
    fxtmp = fytmp = fztmp = 0.0;
    extmp = eytmp = eztmp = 0.0;

    // self term Eq. (55) for I_{ii} and Eq. (52) and in Barros et al
    const double curvature_threshold = sqrt(area[i]);
    if (curvature[i] < curvature_threshold) {
      const double sf = curvature[i] / (4.0 * MY_PIS * curvature_threshold) * area[i] * q[i];
      efield[i][0] = sf * norm[i].x;
      efield[i][1] = sf * norm[i].y;
      efield[i][2] = sf * norm[i].z;
    } else {
      efield[i][0] = efield[i][1] = efield[i][2] = 0;
    }

    epot[i] = 0;

    for (int jj = 0; jj < jnum; jj++) {    // loop over neighbors
      int j = jlist[jj];
      const int ni = sbmask(j);
      const double factor_coul = special_coul[ni];
      j &= NEIGHMASK;

      const double delx = xtmp - x[j].x;
      const double dely = ytmp - x[j].y;
      const double delz = ztmp - x[j].z;
      const double rsq = delx * delx + dely * dely + delz * delz;
      const int jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      const double r2inv = 1.0 / rsq;
      const double r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq)) {    // coulombic
        if (!ncoultablebits || rsq <= tabinnersq) {

          grij = g_ewald * r;
          expm2 = exp(-grij * grij);
          t = 1.0 / (1.0 + EWALD_P * grij);
          erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
          prefactor = qqrd2e * qtmp * q[j] / r;
          force_coul = prefactor * (erfc + EWALD_F * grij * expm2);
          if (factor_coul < 1.0) force_coul -= (1.0 - factor_coul) * prefactor;

          prefactorE = q[j] / r;
          efield_i = prefactorE * (erfc + EWALD_F * grij * expm2);
          if (factor_coul < 1.0) efield_i -= (1.0 - factor_coul) * prefactorE;
          epot_i = efield_i;
        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction * dftable[itable];
          force_coul = qtmp * q[j] * table;
          efield_i = q[j] * table / qqrd2e;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction * dctable[itable];
            prefactor = qtmp * q[j] * table;
            force_coul -= (1.0 - factor_coul) * prefactor;

            prefactorE = q[j] * table / qqrd2e;
            efield_i -= (1.0 - factor_coul) * prefactorE;
          }
          epot_i = efield_i;
        }
      } else
        epot_i = efield_i = force_coul = ecoul = 0.0;

      if (rsq < cut_ljsqi[jtype]) {                          // lj
        if (order6) {                                        // long-range lj
          if (!ndisptablebits || rsq <= tabinnerdispsq) {    // series real space
            double rn = r2inv * r2inv * r2inv;
            double x2 = g2 * rsq, a2 = 1.0 / x2;
            x2 = a2 * exp(-x2) * lj4i[jtype];
            if (ni == 0) {
              force_lj = (rn *= rn) * lj1i[jtype] -
                  g8 * (((6.0 * a2 + 6.0) * a2 + 3.0) * a2 + 1.0) * x2 * rsq;
              if (EFLAG) evdwl = rn * lj3i[jtype] - g6 * ((a2 + 1.0) * a2 + 0.5) * x2;
            } else {    // special case
              const double fs = special_lj[ni], ts = rn * (1.0 - fs);
              force_lj = fs * (rn *= rn) * lj1i[jtype] -
                  g8 * (((6.0 * a2 + 6.0) * a2 + 3.0) * a2 + 1.0) * x2 * rsq + ts * lj2i[jtype];
              if (EFLAG)
                evdwl =
                    fs * rn * lj3i[jtype] - g6 * ((a2 + 1.0) * a2 + 0.5) * x2 + ts * lj4i[jtype];
            }
          } else {    // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask) >> ndispshiftbits;
            const double f_disp = (rsq - rdisptable[disp_k]) * drdisptable[disp_k];
            double rn = r2inv * r2inv * r2inv;
            if (ni == 0) {
              force_lj = (rn *= rn) * lj1i[jtype] -
                  (fdisptable[disp_k] + f_disp * dfdisptable[disp_k]) * lj4i[jtype];
              if (EFLAG)
                evdwl = rn * lj3i[jtype] -
                    (edisptable[disp_k] + f_disp * dedisptable[disp_k]) * lj4i[jtype];
            } else {    // special case
              const double fs = special_lj[ni], ts = rn * (1.0 - fs);
              force_lj = fs * (rn *= rn) * lj1i[jtype] -
                  (fdisptable[disp_k] + f_disp * dfdisptable[disp_k]) * lj4i[jtype] +
                  ts * lj2i[jtype];
              if (EFLAG)
                evdwl = fs * rn * lj3i[jtype] -
                    (edisptable[disp_k] + f_disp * dedisptable[disp_k]) * lj4i[jtype] +
                    ts * lj4i[jtype];
            }
          }
        } else {    // cut lj
          const double rn = r2inv * r2inv * r2inv;
          if (ni == 0) {
            force_lj = rn * (rn * lj1i[jtype] - lj2i[jtype]);
            if (EFLAG) evdwl = rn * (rn * lj3i[jtype] - lj4i[jtype]) - offseti[jtype];
          } else {    // special case
            const double fs = special_lj[ni];
            force_lj = fs * rn * (rn * lj1i[jtype] - lj2i[jtype]);
            if (EFLAG) evdwl = fs * (rn * (rn * lj3i[jtype] - lj4i[jtype]) - offseti[jtype]);
          }
        }
      } else
        force_lj = evdwl = 0.0;

      const double fpair_i = (force_coul * etmp + force_lj) * r2inv;
      fxtmp += delx * fpair_i;
      fytmp += dely * fpair_i;
      fztmp += delz * fpair_i;

      efield_i *= (etmp * r2inv);
      extmp += delx * efield_i;
      eytmp += dely * efield_i;
      eztmp += delz * efield_i;

      epot[i] += epot_i;

      if (EFLAG) {
        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq)
            ecoul = prefactor * 0.5 * (etmp + eps[j]) * erfc;
          else {
            table = etable[itable] + fraction * detable[itable];
            ecoul = qtmp * q[j] * 0.5 * (etmp + eps[j]) * table;
          }
          if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
        } else
          ecoul = 0.0;
      }

      if (EVFLAG) ev_tally_full_thr(this, i, evdwl, ecoul, fpair_i, delx, dely, delz, thr);
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    efield[i][0] += extmp;
    efield[i][1] += eytmp;
    efield[i][2] += eztmp;
  }
}

/* ---------------------------------------------------------------------- */

double PairLJLongCoulLongDielectricOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJLongCoulLongDielectric::memory_usage();

  return bytes;
}
