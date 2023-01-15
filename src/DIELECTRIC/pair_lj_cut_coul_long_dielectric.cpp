/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Trung Nguyen (Northwestern)
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_long_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace EwaldConst;
using MathConst::MY_PIS;

static constexpr double EPSILON = 1.0e-6;

/* ---------------------------------------------------------------------- */

PairLJCutCoulLongDielectric::PairLJCutCoulLongDielectric(LAMMPS *_lmp) : PairLJCutCoulLong(_lmp)
{
  respa_enable = 0;
  cut_respa = nullptr;
  efield = nullptr;
  epot = nullptr;
  nmax = 0;
  no_virial_fdotr_compute = 1;
}

/* ---------------------------------------------------------------------- */

PairLJCutCoulLongDielectric::~PairLJCutCoulLongDielectric()
{
  memory->destroy(efield);
  memory->destroy(epot);
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulLongDielectric::compute(int eflag, int vflag)
{
  int i, ii, j, jj, inum, jnum, itype, jtype, itable;
  double qtmp, etmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul;
  double fpair_i;
  double fraction, table;
  double r, rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
  double grij, expm2, prefactor, t, erfc, prefactorE, efield_i, epot_i;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    memory->destroy(epot);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
    memory->create(epot, nmax, "pair:epot");
  }

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q_scaled;
  double *eps = atom->epsilon;
  double **norm = atom->mu;
  double *curvature = atom->curvature;
  double *area = atom->area;
  int *type = atom->type;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    etmp = eps[i];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // self term Eq. (55) for I_{ii} and Eq. (52) and in Barros et al.

    double curvature_threshold = sqrt(area[i]);
    if (curvature[i] < curvature_threshold) {
      double sf = curvature[i] / (4.0 * MY_PIS * curvature_threshold) * area[i] * q[i];
      efield[i][0] = sf * norm[i][0];
      efield[i][1] = sf * norm[i][1];
      efield[i][2] = sf * norm[i][2];
    } else {
      efield[i][0] = efield[i][1] = efield[i][2] = 0;
    }

    epot[i] = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        r = sqrt(rsq);

        if (rsq < cut_coulsq && rsq > EPSILON) {
          if (!ncoultablebits || rsq <= tabinnersq) {

            grij = g_ewald * r;
            expm2 = exp(-grij * grij);
            t = 1.0 / (1.0 + EWALD_P * grij);
            erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
            prefactor = qqrd2e * qtmp * q[j] / r;
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

            prefactorE = qqrd2e * q[j] / r;
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
            forcecoul = qtmp * q[j] * table;
            efield_i = q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction * dctable[itable];
              prefactor = qtmp * q[j] * table;
              forcecoul -= (1.0 - factor_coul) * prefactor;

              prefactorE = q[j] * table;
              efield_i -= (1.0 - factor_coul) * prefactorE;
            }
            epot_i = efield_i;
          }
        } else
          epot_i = efield_i = forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
        } else
          forcelj = 0.0;

        fpair_i = (forcecoul * etmp + factor_lj * forcelj) * r2inv;
        f[i][0] += delx * fpair_i;
        f[i][1] += dely * fpair_i;
        f[i][2] += delz * fpair_i;

        efield_i *= (etmp * r2inv);
        efield[i][0] += delx * efield_i;
        efield[i][1] += dely * efield_i;
        efield[i][2] += delz * efield_i;

        epot[i] += epot_i;

        if (eflag) {
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

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (evflag) ev_tally_full(i, evdwl, ecoul, fpair_i, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulLongDielectric::init_style()
{
  avec = dynamic_cast<AtomVecDielectric *>(atom->style_match("dielectric"));
  if (!avec) error->all(FLERR, "Pair lj/cut/coul/long/dielectric requires atom style dielectric");

  neighbor->add_request(this, NeighConst::REQ_FULL);

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, cut_respa);
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulLongDielectric::single(int i, int j, int itype, int jtype, double rsq,
                                           double factor_coul, double factor_lj, double &fforce)
{
  double r2inv, r6inv, r, grij, expm2, t, erfc, ei, ej, prefactor;
  double fraction, table, forcecoul, forcelj, phicoul, philj;
  int itable;
  double *q = atom->q_scaled;
  double *eps = atom->epsilon;

  r2inv = 1.0 / rsq;
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq) {
      r = sqrt(rsq);
      grij = g_ewald * r;
      expm2 = exp(-grij * grij);
      t = 1.0 / (1.0 + EWALD_P * grij);
      erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
      prefactor = force->qqrd2e * q[i] * q[j] / r;
      forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
      if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;
    } else {
      union_int_float_t rsq_lookup_single;
      rsq_lookup_single.f = rsq;
      itable = rsq_lookup_single.i & ncoulmask;
      itable >>= ncoulshiftbits;
      fraction = (rsq_lookup_single.f - rtable[itable]) * drtable[itable];
      table = ftable[itable] + fraction * dftable[itable];
      forcecoul = q[i] * q[j] * table;
      if (factor_coul < 1.0) {
        table = ctable[itable] + fraction * dctable[itable];
        prefactor = q[i] * q[j] * table;
        forcecoul -= (1.0 - factor_coul) * prefactor;
      }
    }
  } else
    forcecoul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv * r2inv * r2inv;
    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
  } else
    forcelj = 0.0;

  fforce = (forcecoul * eps[i] + factor_lj * forcelj) * r2inv;

  double eng = 0.0;
  if (eps[i] == 1)
    ei = 0;
  else
    ei = eps[i];
  if (eps[j] == 1)
    ej = 0;
  else
    ej = eps[j];
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq)
      phicoul = prefactor * 0.5 * (ei + ej) * erfc;
    else {
      table = etable[itable] + fraction * detable[itable];
      phicoul = q[i] * q[j] * 0.5 * (ei + ej) * table;
    }
    if (factor_coul < 1.0) phicoul -= (1.0 - factor_coul) * prefactor;
    eng += phicoul;
  }

  if (rsq < cut_ljsq[itype][jtype]) {
    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
    eng += factor_lj * philj;
  }

  return eng;
}
