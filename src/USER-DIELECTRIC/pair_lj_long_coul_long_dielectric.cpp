/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_lj_long_coul_long_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "integrate.h"
#include "kspace.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace MathExtra;

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

/* ---------------------------------------------------------------------- */

PairLJLongCoulLongDielectric::PairLJLongCoulLongDielectric(LAMMPS *lmp) : PairLJLongCoulLong(lmp)
{
  respa_enable = 0;
  cut_respa = nullptr;
  efield = nullptr;
  epot = nullptr;
  nmax = 0;
}

/* ----------------------------------------------------------------------
   free all arrays
------------------------------------------------------------------------- */

PairLJLongCoulLongDielectric::~PairLJLongCoulLongDielectric()
{
  memory->destroy(efield);
  memory->destroy(epot);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJLongCoulLongDielectric::init_style()
{
  PairLJLongCoulLong::init_style();

  avec = (AtomVecDielectric *) atom->style_match("dielectric");
  if (!avec) error->all(FLERR, "Pair lj/long/coul/long/dielectric requires atom style dielectric");

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ----------------------------------------------------------------------
   compute pair interactions
------------------------------------------------------------------------- */

void PairLJLongCoulLongDielectric::compute(int eflag, int vflag)
{
  double evdwl, ecoul;
  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    memory->destroy(epot);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
    memory->create(epot, nmax, "pair:epot");
  }

  double **x = atom->x, *x0 = x[0];
  double **f = atom->f;
  double *q = atom->q;
  double *eps = avec->epsilon;
  double **norm = avec->mu;
  double *curvature = avec->curvature;
  double *area = avec->area;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  int i, j, itype, jtype, itable;
  double qtmp, etmp, xtmp, ytmp, ztmp, delx, dely, delz;
  int order1 = ewald_order & (1 << 1), order6 = ewald_order & (1 << 6);
  int *ineigh, *ineighn, *jneigh, *jneighn, ni;
  double fpair_i, fpair_j;
  double fraction, table;
  double *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double grij, expm2, prefactor, t, erfc, prefactorE, efield_i, epot_i;
  double r, rsq, r2inv, force_coul, force_lj, factor_coul, factor_lj;
  double g2 = g_ewald_6 * g_ewald_6, g6 = g2 * g2 * g2, g8 = g6 * g2;
  double xi[3];

  ineighn = (ineigh = list->ilist) + list->inum;

  for (; ineigh < ineighn; ++ineigh) {    // loop over my atoms
    i = *ineigh;
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    etmp = eps[i];

    offseti = offset[itype = type[i]];
    lj1i = lj1[itype];
    lj2i = lj2[itype];
    lj3i = lj3[itype];
    lj4i = lj4[itype];
    cut_ljsqi = cut_ljsq[itype];
    memcpy(xi, x0 + (i + (i << 1)), 3 * sizeof(double));
    jneighn = (jneigh = list->firstneigh[i]) + list->numneigh[i];

    // self term Eq. (55) for I_{ii} and Eq. (52) and in Barros et al
    double curvature_threshold = sqrt(area[i]);
    if (curvature[i] < curvature_threshold) {
      double sf = curvature[i] / (4.0 * MY_PIS * curvature_threshold) * area[i] * q[i];
      efield[i][0] = sf * norm[i][0];
      efield[i][1] = sf * norm[i][1];
      efield[i][2] = sf * norm[i][2];
    } else {
      efield[i][0] = efield[i][1] = efield[i][2] = 0;
    }

    epot[i] = 0;

    for (; jneigh < jneighn; ++jneigh) {    // loop over neighbors
      j = *jneigh;
      ni = sbmask(j);
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq >= cutsq[itype][jtype]) continue;

      r2inv = 1.0 / rsq;
      r = sqrt(rsq);

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
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
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
              if (eflag) evdwl = rn * lj3i[jtype] - g6 * ((a2 + 1.0) * a2 + 0.5) * x2;
            } else {    // special case
              double f = special_lj[ni], t = rn * (1.0 - f);
              force_lj = f * (rn *= rn) * lj1i[jtype] -
                  g8 * (((6.0 * a2 + 6.0) * a2 + 3.0) * a2 + 1.0) * x2 * rsq + t * lj2i[jtype];
              if (eflag)
                evdwl = f * rn * lj3i[jtype] - g6 * ((a2 + 1.0) * a2 + 0.5) * x2 + t * lj4i[jtype];
            }
          } else {    // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask) >> ndispshiftbits;
            double f_disp = (rsq - rdisptable[disp_k]) * drdisptable[disp_k];
            double rn = r2inv * r2inv * r2inv;
            if (ni == 0) {
              force_lj = (rn *= rn) * lj1i[jtype] -
                  (fdisptable[disp_k] + f_disp * dfdisptable[disp_k]) * lj4i[jtype];
              if (eflag)
                evdwl = rn * lj3i[jtype] -
                    (edisptable[disp_k] + f_disp * dedisptable[disp_k]) * lj4i[jtype];
            } else {    // special case
              double f = special_lj[ni], t = rn * (1.0 - f);
              force_lj = f * (rn *= rn) * lj1i[jtype] -
                  (fdisptable[disp_k] + f_disp * dfdisptable[disp_k]) * lj4i[jtype] +
                  t * lj2i[jtype];
              if (eflag)
                evdwl = f * rn * lj3i[jtype] -
                    (edisptable[disp_k] + f_disp * dedisptable[disp_k]) * lj4i[jtype] +
                    t * lj4i[jtype];
            }
          }
        } else {    // cut lj
          double rn = r2inv * r2inv * r2inv;
          if (ni == 0) {
            force_lj = rn * (rn * lj1i[jtype] - lj2i[jtype]);
            if (eflag) evdwl = rn * (rn * lj3i[jtype] - lj4i[jtype]) - offseti[jtype];
          } else {    // special case
            double f = special_lj[ni];
            force_lj = f * rn * (rn * lj1i[jtype] - lj2i[jtype]);
            if (eflag) evdwl = f * (rn * (rn * lj3i[jtype] - lj4i[jtype]) - offseti[jtype]);
          }
        }
      } else force_lj = evdwl = 0.0;

      fpair_i = (force_coul * etmp + force_lj) * r2inv;
      f[i][0] += delx * fpair_i;
      f[i][1] += dely * fpair_i;
      f[i][2] += delz * fpair_i;

      efield_i *= (etmp * r2inv);
      efield[i][0] += delx * efield_i;
      efield[i][1] += dely * efield_i;
      efield[i][2] += delz * efield_i;

      epot[i] += epot_i;

      if (newton_pair && j >= nlocal) {

        fpair_j = (force_coul * eps[j] + factor_lj * force_lj) * r2inv;
        f[j][0] -= delx * fpair_j;
        f[j][1] -= dely * fpair_j;
        f[j][2] -= delz * fpair_j;
      }

      if (eflag) {
        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq)
            ecoul = prefactor * (etmp + eps[j]) * erfc;
          else {
            table = etable[itable] + fraction * detable[itable];
            ecoul = qtmp * q[j] * (etmp + eps[j]) * table;
          }
          ecoul *= 0.5;
          if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
        } else
          ecoul = 0.0;
      }

      if (evflag) ev_tally_full(i, evdwl, ecoul, fpair_i, delx, dely, delz);
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairLJLongCoulLongDielectric::single(int i, int j, int itype, int jtype, double rsq,
                                            double factor_coul, double factor_lj, double &fforce)
{
  double r2inv, r6inv, force_coul, force_lj, ei, ej;
  double g2 = g_ewald_6 * g_ewald_6, g6 = g2 * g2 * g2, g8 = g6 * g2, *q = atom->q;
  double *eps = avec->epsilon;

  double eng = 0.0;
  if (eps[i] == 1)
    ei = 0;
  else
    ei = eps[i];
  if (eps[j] == 1)
    ej = 0;
  else
    ej = eps[j];

  r2inv = 1.0 / rsq;
  if ((ewald_order & 2) && (rsq < cut_coulsq)) {    // coulombic
    if (!ncoultablebits || rsq <= tabinnersq) {     // series real space
      double r = sqrt(rsq), x = g_ewald * r;
      double s = force->qqrd2e * q[i] * q[j], t = 1.0 / (1.0 + EWALD_P * x);
      r = s * (1.0 - factor_coul) / r;
      s *= g_ewald * exp(-x * x);
      force_coul = (t *= ((((t * A5 + A4) * t + A3) * t + A2) * t + A1) * s / x) + EWALD_F * s - r;
      eng += (t - r) * (ei + ej) * 0.5;
    } else {    // table real space
      union_int_float_t t;
      t.f = rsq;
      const int k = (t.i & ncoulmask) >> ncoulshiftbits;
      double f = (rsq - rtable[k]) * drtable[k], qiqj = q[i] * q[j];
      t.f = (1.0 - factor_coul) * (ctable[k] + f * dctable[k]);
      force_coul = qiqj * (ftable[k] + f * dftable[k] - t.f);
      eng += qiqj * (etable[k] + f * detable[k] - t.f) * (ei + ej) * 0.5;
    }
  } else
    force_coul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {    // lennard-jones
    r6inv = r2inv * r2inv * r2inv;
    if (ewald_order & 64) {    // long-range
      double x2 = g2 * rsq, a2 = 1.0 / x2, t = r6inv * (1.0 - factor_lj);
      x2 = a2 * exp(-x2) * lj4[itype][jtype];
      force_lj = factor_lj * (r6inv *= r6inv) * lj1[itype][jtype] -
          g8 * (((6.0 * a2 + 6.0) * a2 + 3.0) * a2 + a2) * x2 * rsq + t * lj2[itype][jtype];
      eng += factor_lj * r6inv * lj3[itype][jtype] - g6 * ((a2 + 1.0) * a2 + 0.5) * x2 +
          t * lj4[itype][jtype];
    } else {    // cut
      force_lj = factor_lj * r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
      eng += factor_lj *
          (r6inv * (r6inv * lj3[itype][jtype] - lj4[itype][jtype]) - offset[itype][jtype]);
    }
  } else
    force_lj = 0.0;

  fforce = (force_coul * eps[i] + force_lj) * r2inv;
  return eng;
}
