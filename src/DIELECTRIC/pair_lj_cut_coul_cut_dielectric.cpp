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

#include "pair_lj_cut_coul_cut_dielectric.h"

#include "atom.h"
#include "atom_vec_dielectric.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"

#include <cmath>

using namespace LAMMPS_NS;
using namespace MathConst;

#define EPSILON 1e-6

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutDielectric::PairLJCutCoulCutDielectric(LAMMPS *lmp) : PairLJCutCoulCut(lmp)
{
  efield = nullptr;
  epot = nullptr;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutDielectric::~PairLJCutCoulCutDielectric()
{
  memory->destroy(efield);
  memory->destroy(epot);
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutDielectric::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double qtmp, etmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul;
  double fpair_i, fpair_j;
  double rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj, efield_i, epot_i;
  int *ilist, *jlist, *numneigh, **firstneigh;

  if (atom->nmax > nmax) {
    memory->destroy(efield);
    memory->destroy(epot);
    nmax = atom->nmax;
    memory->create(efield, nmax, 3, "pair:efield");
    memory->create(epot, nmax, "pair:epot");
  }

  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double *eps = atom->epsilon;
  double **norm = atom->mu;
  double *curvature = atom->curvature;
  double *area = atom->area;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
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
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    etmp = eps[i];
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

    epot[i] = 0;

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

        if (rsq < cut_coulsq[itype][jtype] && rsq > EPSILON) {
          efield_i = q[j] * sqrt(r2inv);
          forcecoul = qqrd2e * qtmp * efield_i;
          epot_i = efield_i;
        } else
          epot_i = efield_i = forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
        } else
          forcelj = 0.0;

        fpair_i = (factor_coul * etmp * forcecoul + factor_lj * forcelj) * r2inv;
        f[i][0] += delx * fpair_i;
        f[i][1] += dely * fpair_i;
        f[i][2] += delz * fpair_i;

        efield_i *= (factor_coul * etmp * r2inv);
        efield[i][0] += delx * efield_i;
        efield[i][1] += dely * efield_i;
        efield[i][2] += delz * efield_i;

        epot[i] += epot_i;

        if (newton_pair && j >= nlocal) {
          fpair_j = (factor_coul * eps[j] * forcecoul + factor_lj * forcelj) * r2inv;
          f[j][0] -= delx * fpair_j;
          f[j][1] -= dely * fpair_j;
          f[j][2] -= delz * fpair_j;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype]) {
            ecoul = factor_coul * qqrd2e * qtmp * q[j] * (etmp + eps[j]) * sqrt(r2inv);
          } else
            ecoul = 0.0;
          ecoul *= 0.5;
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

void PairLJCutCoulCutDielectric::init_style()
{
  avec = (AtomVecDielectric *) atom->style_match("dielectric");
  if (!avec) error->all(FLERR, "Pair lj/cut/coul/cut/dielectric requires atom style dielectric");

  int irequest = neighbor->request(this, instance_me);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1;
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulCutDielectric::single(int i, int j, int itype, int jtype, double rsq,
                                          double factor_coul, double factor_lj, double &fforce)
{
  double r2inv, r6inv, forcecoul, forcelj, phicoul, ei, ej, philj;
  double *eps = atom->epsilon;

  r2inv = 1.0 / rsq;
  if (rsq < cut_coulsq[itype][jtype])
    forcecoul = force->qqrd2e * atom->q[i] * atom->q[j] * sqrt(r2inv) * eps[i];
  else
    forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv * r2inv * r2inv;
    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
  } else
    forcelj = 0.0;
  fforce = (factor_coul * forcecoul + factor_lj * forcelj) * r2inv;

  double eng = 0.0;
  if (eps[i] == 1)
    ei = 0;
  else
    ei = eps[i];
  if (eps[j] == 1)
    ej = 0;
  else
    ej = eps[j];
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i] * atom->q[j] * sqrt(r2inv);
    phicoul *= 0.5 * (ei + ej);
    eng += factor_coul * phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
    eng += factor_lj * philj;
  }

  return eng;
}
