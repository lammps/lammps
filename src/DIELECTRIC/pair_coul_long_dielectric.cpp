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

#include "pair_coul_long_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

/* ---------------------------------------------------------------------- */

PairCoulLongDielectric::PairCoulLongDielectric(LAMMPS *lmp) : PairCoulLong(lmp)
{
  efield = nullptr;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

PairCoulLongDielectric::~PairCoulLongDielectric()
{
  memory->destroy(efield);
}

/* ---------------------------------------------------------------------- */

void PairCoulLongDielectric::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itable, itype, jtype;
  double qtmp, etmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul;
  double fpair_i, fpair_j;
  double fraction, table;
  double r, r2inv, forcecoul, factor_coul;
  double grij, expm2, prefactor, t, erfc, prefactorE, efield_i;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double rsq;

  ecoul = 0.0;
  ev_init(eflag, vflag);

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
  }

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double **norm = atom->mu;
  double *curvature = atom->curvature;
  double *area = atom->area;
  double *eps = atom->epsilon;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    etmp = eps[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

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

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {
        r2inv = 1.0 / rsq;
        if (!ncoultablebits || rsq <= tabinnersq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij * grij);
          t = 1.0 / (1.0 + EWALD_P * grij);
          erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
          prefactor = qqrd2e * scale[itype][jtype] * qtmp * q[j] / r;
          forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
          if (factor_coul < 1.0) forcecoul -= (1.0 - factor_coul) * prefactor;

          prefactorE = qqrd2e * scale[itype][jtype] * q[j] / r;
          efield_i = prefactorE * (erfc + EWALD_F * grij * expm2);
          if (factor_coul < 1.0) efield_i -= (1.0 - factor_coul) * prefactorE;

        } else {
          union_int_float_t rsq_lookup;
          rsq_lookup.f = rsq;
          itable = rsq_lookup.i & ncoulmask;
          itable >>= ncoulshiftbits;
          fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
          table = ftable[itable] + fraction * dftable[itable];
          forcecoul = scale[itype][jtype] * qtmp * q[j] * table;
          efield_i = scale[itype][jtype] * q[j] * table;
          if (factor_coul < 1.0) {
            table = ctable[itable] + fraction * dctable[itable];
            prefactor = scale[itype][jtype] * qtmp * q[j] * table;
            forcecoul -= (1.0 - factor_coul) * prefactor;

            prefactorE = scale[itype][jtype] * q[j] * table;
            efield_i -= (1.0 - factor_coul) * prefactorE;
          }
        }

        fpair_i = etmp * forcecoul * r2inv;
        f[i][0] += delx * fpair_i;
        f[i][1] += dely * fpair_i;
        f[i][2] += delz * fpair_i;

        efield_i *= (etmp * r2inv);
        efield[i][0] += delx * efield_i;
        efield[i][1] += dely * efield_i;
        efield[i][2] += delz * efield_i;

        if (newton_pair && j >= nlocal) {
          fpair_j = eps[j] * forcecoul * r2inv;
          f[j][0] -= delx * fpair_j;
          f[j][1] -= dely * fpair_j;
          f[j][2] -= delz * fpair_j;
        }

        if (eflag) {
          if (!ncoultablebits || rsq <= tabinnersq)
            ecoul = prefactor * (etmp + eps[j]) * erfc;
          else {
            table = etable[itable] + fraction * detable[itable];
            ecoul = scale[itype][jtype] * qtmp * q[j] * (etmp + eps[j]) * table;
          }
          ecoul *= 0.5;
          if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor;
        }

        if (evflag) ev_tally_full(i, 0.0, ecoul, fpair_i, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulLongDielectric::init_style()
{
  avec = (AtomVecDielectric *) atom->style_match("dielectric");
  if (!avec) error->all(FLERR, "Pair coul/long/dielectric requires atom style dielectric");

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables

  if (ncoultablebits) init_tables(cut_coul, nullptr);
}
