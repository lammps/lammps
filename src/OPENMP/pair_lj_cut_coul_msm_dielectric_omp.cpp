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

#include "pair_lj_cut_coul_msm_dielectric_omp.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>

#include "omp_compat.h"

using namespace LAMMPS_NS;
using MathConst::MY_PIS;

static constexpr double EPSILON = 1.0e-6;

namespace {
using dbl3_t = struct {
  double x, y, z;
};
}    // namespace

/* ---------------------------------------------------------------------- */

PairLJCutCoulMSMDielectricOMP::PairLJCutCoulMSMDielectricOMP(LAMMPS *_lmp) :
    PairLJCutCoulMSMDielectric(_lmp), ThrOMP(_lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulMSMDielectricOMP::compute(int eflag, int vflag)
{
  if (force->kspace->scalar_pressure_flag)
    error->all(FLERR, "Must use 'kspace_modify pressure/scalar no' with OMP MSM Pair styles");

  ev_init(eflag, vflag);

  if (!efield || atom->nmax > nmax) {
    if (efield) memory->destroy(efield);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
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
void PairLJCutCoulMSMDielectricOMP::eval(int iifrom, int iito, ThrData *const thr)
{
  int i, j, ii, jj, jnum, itype, jtype, itable;
  double qtmp, etmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair_i;
  double fraction, table;
  double r, rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
  double egamma, fgamma, prefactor, prefactorE, efield_i;

  evdwl = ecoul = 0.0;

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
  const double qqrd2e = force->qqrd2e;
  double fxtmp, fytmp, fztmp, extmp, eytmp, eztmp;

  const int *ilist = list->ilist;
  const int *numneigh = list->numneigh;
  int *const *firstneigh = list->firstneigh;

  for (ii = iifrom; ii < iito; ++ii) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i].x;
    ytmp = x[i].y;
    ztmp = x[i].z;
    etmp = eps[i];
    itype = type[i];
    const double *_noalias const cutsqi    = cutsq[itype];
    const double *_noalias const cut_ljsqi = cut_ljsq[itype];
    const double *_noalias const lj1i      = lj1[itype];
    const double *_noalias const lj2i      = lj2[itype];
    const double *_noalias const lj3i      = lj3[itype];
    const double *_noalias const lj4i      = lj4[itype];
    const double *_noalias const offseti   = offset[itype];
    const int *_noalias const jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;
    extmp = eytmp = eztmp = 0.0;

    // self term Eq. (55) for I_{ii} and Eq. (52) in Barros et al.
    double curvature_threshold = sqrt(area[i]);
    if (curvature[i] < curvature_threshold) {
      double sf = curvature[i] / (4.0 * MY_PIS * curvature_threshold) * area[i] * q[i];
      efield[i][0] = sf * norm[i].x;
      efield[i][1] = sf * norm[i].y;
      efield[i][2] = sf * norm[i].z;
    } else {
      efield[i][0] = efield[i][1] = efield[i][2] = 0;
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j].x;
      dely = ytmp - x[j].y;
      delz = ztmp - x[j].z;
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsqi[jtype]) {
        r2inv = 1.0 / rsq;
        prefactor = prefactorE = egamma = r6inv = 0.0;

        if (rsq < cut_coulsq && rsq > EPSILON) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            prefactor = qqrd2e * qtmp * q[j] / r;
            egamma = 1.0 - (r / cut_coul) * force->kspace->gamma(r / cut_coul);
            fgamma = 1.0 + (rsq / cut_coulsq) * force->kspace->dgamma(r / cut_coul);
            forcecoul = prefactor * fgamma;
            if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

            prefactorE = qqrd2e * q[j] / r;
            efield_i = prefactorE * fgamma;
            if (factor_coul < 1.0) efield_i -= (1.0 - factor_coul) * prefactorE;
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = ((double) rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction * dftable[itable];
            forcecoul = qtmp * q[j] * table;
            efield_i = q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction * dctable[itable];
              prefactor = qtmp * q[j] * table;
              forcecoul -= (1.0 - factor_coul) * prefactor;
              prefactorE = q[j] * table;
              efield_i -= (1.0 - factor_coul) * prefactorE;
            }
          }
        } else
          forcecoul = efield_i = 0.0;

        if (rsq < cut_ljsqi[jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1i[jtype] * r6inv - lj2i[jtype]);
        } else
          forcelj = 0.0;

        fpair_i = (forcecoul * etmp + factor_lj * forcelj) * r2inv;
        fxtmp += delx * fpair_i;
        fytmp += dely * fpair_i;
        fztmp += delz * fpair_i;

        efield_i *= (etmp * r2inv);
        extmp += delx * efield_i;
        eytmp += dely * efield_i;
        eztmp += delz * efield_i;

        if (EFLAG) {
          if (rsq < cut_coulsq && rsq > EPSILON) {
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor * 0.5 * (etmp + eps[j]) * egamma;
            else {
              table = etable[itable] + fraction * detable[itable];
              ecoul = qtmp * q[j] * 0.5 * (etmp + eps[j]) * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
          } else
            ecoul = 0.0;

          if (rsq < cut_ljsqi[jtype]) {
            evdwl = r6inv * (lj3i[jtype] * r6inv - lj4i[jtype]) - offseti[jtype];
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (EVFLAG) ev_tally_full_thr(this, i, evdwl, ecoul, fpair_i, delx, dely, delz, thr);
      }
    }
    f[i].x += fxtmp;
    f[i].y += fytmp;
    f[i].z += fztmp;
    efield[i][0] += extmp;
    efield[i][1] += eytmp;
    efield[i][2] += eztmp;
  }
}
