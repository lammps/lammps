/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mateo Rodríguez (mateorsuarez@gmail.com) (IFF-CSIC)
   Work done at the Molecular Interactions Group (INTERMOL) of the
   Fundamental Physics Institute (http://intermol.iff.csic.es/).
   Optimization of the code: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_lj_pirani.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathSpecial::square;

/* ---------------------------------------------------------------------- */

PairLJPirani::PairLJPirani(LAMMPS *lmp) : Pair(lmp), cut_respa(nullptr)
{
  respa_enable = 1;
  born_matrix_enable = 0;
  writedata = 1;
}
/* ---------------------------------------------------------------------- */

PairLJPirani::~PairLJPirani()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(alpha);
    memory->destroy(beta);
    memory->destroy(gamma);
    memory->destroy(rm);
    memory->destroy(epsilon);
    memory->destroy(offset);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJPirani::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n, n, "pair:cutsq");

  memory->create(cut, n, n, "pair:cut");
  memory->create(alpha, n, n, "pair:alpha");
  memory->create(beta, n, n, "pair:beta");
  memory->create(gamma, n, n, "pair:gamma");
  memory->create(rm, n, n, "pair:rm");
  memory->create(epsilon, n, n, "pair:epsilon");
  memory->create(offset, n, n, "pair:offset");
}

/* ---------------------------------------------------------------------- */

void PairLJPirani::compute(int eflag, int vflag)

{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double r, rx, n_x;
  double pow_rx_n_x, pow_rx_gamma;
  double filj1, filj2, filj3, filj4, filj5, filj6, forceilj;
  double ilj1, ilj2;
  double fxtmp, fytmp, fztmp;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        r = sqrt(rsq);

        rx = r / rm[itype][jtype];
        n_x = alpha[itype][jtype] * rx * rx + beta[itype][jtype];
        pow_rx_n_x = pow(1.0 / rx, n_x);
        pow_rx_gamma = pow(1.0 / rx, gamma[itype][jtype]);

        filj1 = -2.0 * alpha[itype][jtype] * gamma[itype][jtype] * rx * pow_rx_n_x /
            (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

        filj2 = +2.0 * alpha[itype][jtype] * rx * n_x * pow_rx_gamma /
            (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

        filj3 = -2.0 * alpha[itype][jtype] * rx * pow_rx_gamma /
            (rm[itype][jtype] * (n_x - gamma[itype][jtype]));

        filj4 = +2.0 * alpha[itype][jtype] * gamma[itype][jtype] * (rx / rm[itype][jtype]) *
            log(1 / rx) * pow_rx_n_x / (n_x - gamma[itype][jtype]);

        filj5 = -1.0 * gamma[itype][jtype] * n_x * pow_rx_n_x / (r * (n_x - gamma[itype][jtype]));

        filj6 = +1.0 * gamma[itype][jtype] * n_x * pow_rx_gamma / (r * (n_x - gamma[itype][jtype]));

        // F = -dV/dr
        forceilj = -epsilon[itype][jtype] * (filj1 + filj2 + filj3 + filj4 + filj5 + filj6);
        fpair = factor_lj * forceilj / r;    // F_x = -x/r * dV/dr (chain rule)

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          ilj1 = epsilon[itype][jtype] * gamma[itype][jtype] * pow(1 / rx, n_x) /
              (n_x - gamma[itype][jtype]);
          ilj2 = -epsilon[itype][jtype] * n_x * pow(1 / rx, gamma[itype][jtype]) /
              (n_x - gamma[itype][jtype]);

          evdwl = ilj1 + ilj2 - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void PairLJPirani::compute_inner()
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double rsq, factor_lj, rsw;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double r, rx, n_x;
  double pow_rx_n_x, pow_rx_gamma;
  double filj1, filj2, filj3, filj4, filj5, filj6, forceilj;
  double fxtmp, fytmp, fztmp;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum_inner;
  ilist = list->ilist_inner;
  numneigh = list->numneigh_inner;
  firstneigh = list->firstneigh_inner;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];

  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on * cut_out_on;
  double cut_out_off_sq = cut_out_off * cut_out_off;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_out_off_sq) {
        jtype = type[j];
        r = sqrt(rsq);

        rx = r / rm[itype][jtype];
        n_x = alpha[itype][jtype] * rx * rx + beta[itype][jtype];
        pow_rx_n_x = pow(1.0 / rx, n_x);
        pow_rx_gamma = pow(1.0 / rx, gamma[itype][jtype]);

        filj1 = -2.0 * alpha[itype][jtype] * gamma[itype][jtype] * rx * pow_rx_n_x /
            (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

        filj2 = +2.0 * alpha[itype][jtype] * rx * n_x * pow_rx_gamma /
            (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

        filj3 = -2.0 * alpha[itype][jtype] * rx * pow_rx_gamma /
            (rm[itype][jtype] * (n_x - gamma[itype][jtype]));

        filj4 = +2.0 * alpha[itype][jtype] * gamma[itype][jtype] * (rx / rm[itype][jtype]) *
            log(1 / rx) * pow_rx_n_x / (n_x - gamma[itype][jtype]);

        filj5 = -1.0 * gamma[itype][jtype] * n_x * pow_rx_n_x / (r * (n_x - gamma[itype][jtype]));

        filj6 = +1.0 * gamma[itype][jtype] * n_x * pow_rx_gamma / (r * (n_x - gamma[itype][jtype]));

        // F = -dV/dr
        forceilj = -epsilon[itype][jtype] * (filj1 + filj2 + filj3 + filj4 + filj5 + filj6);
        fpair = factor_lj * forceilj / r;    // F_x = -x/r * dV/dr (chain rule)

        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on) / cut_out_diff;
          fpair *= 1.0 - rsw * rsw * (3.0 - 2.0 * rsw);
        }

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJPirani::compute_middle()

{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, fpair;
  double rsq, factor_lj, rsw;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double r, rx, n_x;
  double pow_rx_n_x, pow_rx_gamma;
  double filj1, filj2, filj3, filj4, filj5, filj6, forceilj;
  double fxtmp, fytmp, fztmp;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum_middle;
  ilist = list->ilist_middle;
  numneigh = list->numneigh_middle;
  firstneigh = list->firstneigh_middle;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off * cut_in_off;
  double cut_in_on_sq = cut_in_on * cut_in_on;
  double cut_out_on_sq = cut_out_on * cut_out_on;
  double cut_out_off_sq = cut_out_off * cut_out_off;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
        jtype = type[j];
        r = sqrt(rsq);

        rx = r / rm[itype][jtype];
        n_x = alpha[itype][jtype] * rx * rx + beta[itype][jtype];
        pow_rx_n_x = pow(1.0 / rx, n_x);
        pow_rx_gamma = pow(1.0 / rx, gamma[itype][jtype]);

        filj1 = -2.0 * alpha[itype][jtype] * gamma[itype][jtype] * rx * pow_rx_n_x /
            (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

        filj2 = +2.0 * alpha[itype][jtype] * rx * n_x * pow_rx_gamma /
            (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

        filj3 = -2.0 * alpha[itype][jtype] * rx * pow_rx_gamma /
            (rm[itype][jtype] * (n_x - gamma[itype][jtype]));

        filj4 = +2.0 * alpha[itype][jtype] * gamma[itype][jtype] * (rx / rm[itype][jtype]) *
            log(1 / rx) * pow_rx_n_x / (n_x - gamma[itype][jtype]);

        filj5 = -1.0 * gamma[itype][jtype] * n_x * pow_rx_n_x / (r * (n_x - gamma[itype][jtype]));

        filj6 = +1.0 * gamma[itype][jtype] * n_x * pow_rx_gamma / (r * (n_x - gamma[itype][jtype]));

        // F = -dV/dr
        forceilj = -epsilon[itype][jtype] * (filj1 + filj2 + filj3 + filj4 + filj5 + filj6);
        fpair = factor_lj * forceilj / r;
        if (rsq < cut_in_on_sq) {
          rsw = (sqrt(rsq) - cut_in_off) / cut_in_diff;
          fpair *= rsw * rsw * (3.0 - 2.0 * rsw);
        }
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on) / cut_out_diff;
          fpair *= 1.0 + rsw * rsw * (2.0 * rsw - 3.0);
        }

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJPirani::compute_outer(int eflag, int vflag)

{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, factor_lj, rsw;
  int *ilist, *jlist, *numneigh, **firstneigh;

  double r, rx, n_x;
  double pow_rx_n_x, pow_rx_gamma;
  double filj1, filj2, filj3, filj4, filj5, filj6, forceilj;
  double ilj1, ilj2;
  double fxtmp, fytmp, fztmp;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off * cut_in_off;
  double cut_in_on_sq = cut_in_on * cut_in_on;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    fxtmp = fytmp = fztmp = 0.0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        if (rsq > cut_in_off_sq) {
          r = sqrt(rsq);

          rx = r / rm[itype][jtype];
          n_x = alpha[itype][jtype] * rx * rx + beta[itype][jtype];
          pow_rx_n_x = pow(1.0 / rx, n_x);
          pow_rx_gamma = pow(1.0 / rx, gamma[itype][jtype]);

          filj1 = -2.0 * alpha[itype][jtype] * gamma[itype][jtype] * rx * pow_rx_n_x /
              (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

          filj2 = +2.0 * alpha[itype][jtype] * rx * n_x * pow_rx_gamma /
              (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

          filj3 = -2.0 * alpha[itype][jtype] * rx * pow_rx_gamma /
              (rm[itype][jtype] * (n_x - gamma[itype][jtype]));

          filj4 = +2.0 * alpha[itype][jtype] * gamma[itype][jtype] * (rx / rm[itype][jtype]) *
              log(1 / rx) * pow_rx_n_x / (n_x - gamma[itype][jtype]);

          filj5 = -1.0 * gamma[itype][jtype] * n_x * pow_rx_n_x / (r * (n_x - gamma[itype][jtype]));

          filj6 =
              +1.0 * gamma[itype][jtype] * n_x * pow_rx_gamma / (r * (n_x - gamma[itype][jtype]));

          // F = -dV/dr
          forceilj = -epsilon[itype][jtype] * (filj1 + filj2 + filj3 + filj4 + filj5 + filj6);
          fpair = factor_lj * forceilj / r;
          if (rsq < cut_in_on_sq) {
            rsw = (sqrt(rsq) - cut_in_off) / cut_in_diff;
            fpair *= rsw * rsw * (3.0 - 2.0 * rsw);
          }

          fxtmp += delx * fpair;
          fytmp += dely * fpair;
          fztmp += delz * fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx * fpair;
            f[j][1] -= dely * fpair;
            f[j][2] -= delz * fpair;
          }
        }

        if (eflag) {

          r = sqrt(rsq);

          rx = r / rm[itype][jtype];
          n_x = alpha[itype][jtype] * rx * rx + beta[itype][jtype];

          ilj1 = epsilon[itype][jtype] * gamma[itype][jtype] * pow(1 / rx, n_x) /
              (n_x - gamma[itype][jtype]);
          ilj2 = -epsilon[itype][jtype] * n_x * pow(1 / rx, gamma[itype][jtype]) /
              (n_x - gamma[itype][jtype]);

          evdwl = ilj1 + ilj2 - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (vflag) {
          if (rsq <= cut_in_off_sq) {

            r = sqrt(rsq);

            rx = r / rm[itype][jtype];
            n_x = alpha[itype][jtype] * rx * rx + beta[itype][jtype];
            pow_rx_n_x = pow(1.0 / rx, n_x);
            pow_rx_gamma = pow(1.0 / rx, gamma[itype][jtype]);

            filj1 = -2.0 * alpha[itype][jtype] * gamma[itype][jtype] * rx * pow_rx_n_x /
                (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

            filj2 = +2.0 * alpha[itype][jtype] * rx * n_x * pow_rx_gamma /
                (square(n_x - gamma[itype][jtype]) * rm[itype][jtype]);

            filj3 = -2.0 * alpha[itype][jtype] * rx * pow_rx_gamma /
                (rm[itype][jtype] * (n_x - gamma[itype][jtype]));

            filj4 = +2.0 * alpha[itype][jtype] * gamma[itype][jtype] * (rx / rm[itype][jtype]) *
                log(1 / rx) * pow_rx_n_x / (n_x - gamma[itype][jtype]);

            filj5 =
                -1.0 * gamma[itype][jtype] * n_x * pow_rx_n_x / (r * (n_x - gamma[itype][jtype]));

            filj6 =
                +1.0 * gamma[itype][jtype] * n_x * pow_rx_gamma / (r * (n_x - gamma[itype][jtype]));

            // F = -dV/dr
            forceilj = -epsilon[itype][jtype] * (filj1 + filj2 + filj3 + filj4 + filj5 + filj6);
            fpair = factor_lj * forceilj / r;

          } else if (rsq < cut_in_on_sq)
            fpair = factor_lj * forceilj / r;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJPirani::settings(int narg, char **arg)
{
  if (narg != 1)
    error->all(FLERR, "Pair style ilj/cut must have exactly one argument: cutoff distance");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
/*

7 or 8 coefficients: 5 for the ILJ, 2 for the pair, 1 for the cutoff (optional)

*/
void PairLJPirani::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double alpha_one = utils::numeric(FLERR, arg[2], false, lmp);
  double beta_one = utils::numeric(FLERR, arg[3], false, lmp);
  double gamma_one = utils::numeric(FLERR, arg[4], false, lmp);
  double rm_one = utils::numeric(FLERR, arg[5], false, lmp);
  double epsilon_one = utils::numeric(FLERR, arg[6], false, lmp);

  double cut_one = cut_global;
  if (narg == 8) cut_one = utils::numeric(FLERR, arg[7], false, lmp);

  if (rm_one <= 0.0 || epsilon_one < 0.0 || gamma_one <= 0.0)
    error->all(FLERR, "Illegal ILJ coefficients");

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      alpha[i][j] = alpha_one;
      beta[i][j] = beta_one;
      gamma[i][j] = gamma_one;
      rm[i][j] = rm_one;
      epsilon[i][j] = epsilon_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  // Initialize symmetric entries
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      alpha[j][i] = alpha[i][j];
      beta[j][i] = beta[i][j];
      gamma[j][i] = gamma[i][j];
      rm[j][i] = rm[i][j];
      epsilon[j][i] = epsilon[i][j];
      cut[j][i] = cut[i][j];
      setflag[j][i] = setflag[i][j];
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJPirani::init_style()
{
  // request regular or rRESPA neighbor list

  int list_style = NeighConst::REQ_DEFAULT;

  if (update->whichflag == 1 && utils::strmatch(update->integrate_style, "^respa")) {
    auto *respa = dynamic_cast<Respa *>(update->integrate);
    if (respa->level_inner >= 0) list_style = NeighConst::REQ_RESPA_INOUT;
    if (respa->level_middle >= 0) list_style = NeighConst::REQ_RESPA_ALL;
  }
  neighbor->add_request(this, list_style);

  // set rRESPA cutoffs

  if (utils::strmatch(update->integrate_style, "^respa") &&
      (dynamic_cast<Respa *>(update->integrate))->level_inner >= 0)
    cut_respa = (dynamic_cast<Respa *>(update->integrate))->cutoff;
  else
    cut_respa = nullptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJPirani::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  if (offset_flag && (cut[i][j] > 0.0)) {
    double r = cut[i][j] / rm[i][j];
    double nx = alpha[i][j] * r * r + beta[i][j];
    offset[i][j] = epsilon[i][j] *
        ((gamma[i][j] / (nx - gamma[i][j])) * pow(1 / r, nx) -
         (nx / (nx - gamma[i][j])) * pow(1 / r, gamma[i][j]));
  } else
    offset[i][j] = 0.0;

  alpha[j][i] = alpha[i][j];
  beta[j][i] = beta[i][j];
  gamma[j][i] = gamma[i][j];
  rm[j][i] = rm[i][j];
  epsilon[j][i] = epsilon[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR, "Pair cutoff < Respa interior cutoff");

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJPirani::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&alpha[i][j], sizeof(double), 1, fp);
        fwrite(&beta[i][j], sizeof(double), 1, fp);
        fwrite(&gamma[i][j], sizeof(double), 1, fp);
        fwrite(&rm[i][j], sizeof(double), 1, fp);
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJPirani::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJPirani::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &alpha[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &beta[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &gamma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &rm[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&alpha[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&beta[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&gamma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&rm[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJPirani::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJPirani::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %g %g\n", i, alpha[i][i], beta[i][i], gamma[i][i], rm[i][i],
            epsilon[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJPirani::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g\n", i, j, alpha[i][j], beta[i][j], gamma[i][j], rm[i][j],
              epsilon[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJPirani::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                            double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r, rx, n_x, filj1, filj2, filj3, filj4, filj5, filj6, forceilj;
  double ilj1, ilj2;

  r = sqrt(rsq);
  rx = r / rm[itype][jtype];
  n_x = alpha[itype][jtype] * rx * rx + beta[itype][jtype];
  filj1 = -2.0 * alpha[itype][jtype] * gamma[itype][jtype] * rx * pow(1 / rx, n_x) /
      (pow(n_x - gamma[itype][jtype], 2.0) * rm[itype][jtype]);

  filj2 = +2.0 * alpha[itype][jtype] * rx * n_x * pow(1 / rx, gamma[itype][jtype]) /
      (pow(n_x - gamma[itype][jtype], 2.0) * rm[itype][jtype]);

  filj3 = -2.0 * alpha[itype][jtype] * rx * pow(1 / rx, gamma[itype][jtype]) /
      (rm[itype][jtype] * (n_x - gamma[itype][jtype]));

  filj4 = +2.0 * alpha[itype][jtype] * gamma[itype][jtype] * (rx / rm[itype][jtype]) * log(1 / rx) *
      pow(1 / rx, n_x) / (n_x - gamma[itype][jtype]);

  filj5 = -1.0 * gamma[itype][jtype] * n_x * pow(1 / rx, n_x) / (r * (n_x - gamma[itype][jtype]));

  filj6 = +1.0 * gamma[itype][jtype] * n_x * pow(1 / rx, gamma[itype][jtype]) /
      (r * (n_x - gamma[itype][jtype]));

  forceilj = -epsilon[itype][jtype] * (filj1 + filj2 + filj3 + filj4 + filj5 + filj6);

  fforce = factor_lj * forceilj / r;

  ilj1 =
      epsilon[itype][jtype] * gamma[itype][jtype] * pow(1 / rx, n_x) / (n_x - gamma[itype][jtype]);
  ilj2 =
      -epsilon[itype][jtype] * n_x * pow(1 / rx, gamma[itype][jtype]) / (n_x - gamma[itype][jtype]);
  return factor_lj * (ilj1 + ilj2 - offset[itype][jtype]);
}

/* ---------------------------------------------------------------------- */

void *PairLJPirani::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "alpha") == 0) return (void *) alpha;
  if (strcmp(str, "beta") == 0) return (void *) beta;
  if (strcmp(str, "gamma") == 0) return (void *) gamma;
  if (strcmp(str, "rm") == 0) return (void *) rm;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  return nullptr;
}
