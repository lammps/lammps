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
   Contributing authors: Ludwig Ahrens-Iwers (TUHH), Shern Tee (GU), Kamila Savvidi (TUHH), Robert Meissner (Hereon, TUHH)
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_wolf_gauss.h"

#include "atom.h"
#include "comm.h"
#include "electrode_math.h"
#include "error.h"
#include "ewald_const.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cassert>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;
using namespace EwaldConst;

/* ---------------------------------------------------------------------- */

PairLJCutCoulWolfGauss::PairLJCutCoulWolfGauss(LAMMPS *lmp) : Pair(lmp)
{
  ncoultablebits = 0;
  single_enable = 0;
  writedata = 1;
  already_warned = false;
}

/* ---------------------------------------------------------------------- */

PairLJCutCoulWolfGauss::~PairLJCutCoulWolfGauss()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(ispoint);
    memory->destroy(eta);
    memory->destroy(eshift_eta);
    memory->destroy(fshift_eta);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulWolfGauss::compl_error_func(double arg, double expm2)
{
  // expm2 = exp(-arg*arg)
  double const t = 1.0 / (1.0 + EWALD_P * arg);
  return t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::compute(int eflag, int vflag)
{
  int i, ii, j, jj, inum, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double r, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
  double grij, expm2, prefactor, erfc;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
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

  // self and shifted Coulombic energy
  double const alpha_cut = alpha * cut_coul;
  double const expm2_cut = exp(-alpha_cut * alpha_cut);
  double const e_shift = compl_error_func(alpha_cut, expm2_cut) / cut_coul;
  double const pre_self_wolf = (e_shift / 2.0 + alpha / MY_PIS) * qqrd2e;
  double f_shift = -(e_shift + 2.0 * alpha / MY_PIS * expm2_cut) / cut_coul;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    bool ipoint = !!ispoint[itype];

    if (eflag) {
      double const q2 = q[i] * q[i];
      double e = -pre_self_wolf * q2;
      if (!ipoint) e += (eshift_eta[itype][itype] / 2.0 + eta[itype][itype] / MY_PIS) * qqrd2e * q2;
      ev_tally(i, i, nlocal, newton_pair, 0., e, 0., 0., 0., 0.);
    }

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
      bool gausscorr = !ipoint || !ispoint[jtype];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        double erfc_eta = 0.;
        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          grij = alpha * r;
          expm2 = exp(-grij * grij);
          erfc = compl_error_func(grij, expm2);
          prefactor = qqrd2e * qtmp * q[j] / r;
          forcecoul = prefactor * (erfc + EWALD_F * grij * expm2 + f_shift * rsq);
          double forcecorr = 0.0;
          if (gausscorr) {
            double etarij = eta[itype][jtype] * r;
            double expm2_eta = exp(-etarij * etarij);
            erfc_eta = compl_error_func(etarij, expm2_eta);
            forcecorr = erfc_eta + EWALD_F * etarij * expm2_eta;
            forcecoul -= prefactor * (forcecorr + fshift_eta[itype][jtype] * rsq);
          }
          if (factor_coul < 1.0) {
            forcecoul -= (1.0 - factor_coul) * prefactor * (1.0 - forcecorr);
            if (gausscorr) {}
          }
        } else
          forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv * r2inv * r2inv;
          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
        } else
          forcelj = 0.0;

        fpair = (forcecoul + factor_lj * forcelj) * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = prefactor * (erfc - e_shift * r);
            if (gausscorr) ecoul -= prefactor * (erfc_eta - eshift_eta[itype][jtype] * r);
            if (factor_coul < 1.0) ecoul -= (1.0 - factor_coul) * prefactor * (1.0 - erfc_eta);
          } else
            ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut_lj, n + 1, n + 1, "pair:cut_lj");
  memory->create(cut_ljsq, n + 1, n + 1, "pair:cut_ljsq");
  memory->create(epsilon, n + 1, n + 1, "pair:epsilon");
  memory->create(sigma, n + 1, n + 1, "pair:sigma");
  memory->create(ispoint, n + 1, "pair:ispoint");
  memory->create(eta, n + 1, n + 1, "pair:eta");
  memory->create(eshift_eta, n + 1, n + 1, "pair:eshift_eta");
  memory->create(fshift_eta, n + 1, n + 1, "pair:fshift_eta");
  memory->create(lj1, n + 1, n + 1, "pair:lj1");
  memory->create(lj2, n + 1, n + 1, "pair:lj2");
  memory->create(lj3, n + 1, n + 1, "pair:lj3");
  memory->create(lj4, n + 1, n + 1, "pair:lj4");
  memory->create(offset, n + 1, n + 1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR, "Illegal pair_style command");

  alpha = utils::numeric(FLERR, arg[0], false, lmp);
  cut_lj_global = utils::numeric(FLERR, arg[1], false, lmp);
  if (narg == 2)
    cut_coul = cut_lj_global;
  else
    cut_coul = utils::numeric(FLERR, arg[2], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  int ispoint_one = !!(strcmp(arg[4], "NULL") == 0);
  double eta_one = 0.;
  if (!ispoint_one) {
    eta_one = utils::numeric(FLERR, arg[4], false, lmp);
    if (eta_one < MY_SQRT2 * alpha)
      error->all(FLERR, "Reciprocal width is too small for damping parameter {}",
                 alpha);    // sign of Coulomb interaction would flip
    if (ilo != jlo || ihi != jhi)
      if (comm->me == 0)
        error->warning(FLERR, "Gaussian parameter eta cannot be set for mixed interactions");
  }
  double cut_lj_one = cut_lj_global;
  if (narg == 6) cut_lj_one = utils::numeric(FLERR, arg[5], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_lj[i][j] = cut_lj_one;
      if (i == j) {
        ispoint[i] = ispoint_one;
        eta[i][i] = eta_one * MY_ISQRT2;
      }
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR, "Pair style lj/cut/coul/wolf/gauss requires atom attribute q");

  // request regular neighbor list
  int list_style = NeighConst::REQ_DEFAULT;
  neighbor->add_request(this, list_style);

  cut_coulsq = cut_coul * cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutCoulWolfGauss::init_one(int i, int j)
{
  if (i != j) {
    bool const ipoint = !!ispoint[i];
    double const ieta = eta[i][i] * MY_SQRT2;
    bool const jpoint = !!ispoint[j];
    double const jeta = eta[j][j] * MY_SQRT2;

    double tmp = 0;
    if (!ipoint && !jpoint) {
      tmp = ieta * jeta / sqrt(ieta * ieta + jeta * jeta);
    } else if (!ipoint) {
      tmp = ieta;
    } else if (!jpoint) {
      tmp = jeta;
    }
    eta[i][j] = tmp;
    eta[j][i] = tmp;
  }
  double etaij = eta[i][j];
  double const eta_cut_coul = etaij * cut_coul;
  double const expm2 = exp(-eta_cut_coul * eta_cut_coul);
  eshift_eta[i][j] = compl_error_func(eta_cut_coul, expm2) / cut_coul;
  eshift_eta[j][i] = eshift_eta[i][j];
  fshift_eta[i][j] = -(eshift_eta[i][j] + 2.0 * etaij / MY_PIS * expm2) / cut_coul;
  fshift_eta[j][i] = fshift_eta[i][j];

  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut_lj[i][j] = mix_distance(cut_lj[i][i], cut_lj[j][j]);
  }

  double cut = MAX(cut_lj[i][j], cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2], all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

    double sig2 = sigma[i][j] * sigma[i][j];
    double sig6 = sig2 * sig2 * sig2;
    double rc3 = cut_lj[i][j] * cut_lj[i][j] * cut_lj[i][j];
    double rc6 = rc3 * rc3;
    double rc9 = rc3 * rc6;
    etail_ij =
        8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 * (sig6 - 3.0 * rc6) / (9.0 * rc9);
    ptail_ij = 16.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 * (2.0 * sig6 - 3.0 * rc6) /
        (9.0 * rc9);
  }

  return cut;
}

/* ----------------------------------------------------------------------
   compute pair interaction term of vector in constant potential method
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::compute_vector(double *vec, int groupbit, int source_grpbit, bool inv)
{
  point_in_sensor_warning(groupbit);
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int *mask = atom->mask;
  double *special_coul = force->special_coul;
  int const nlocal = atom->nlocal;
  int const inum = list->inum;
  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  double const e_shift = ElectrodeMath::safe_erfc(alpha * cut_coul) / cut_coul;

  for (int ii = 0; ii < inum; ii++) {
    int const i = ilist[ii];
    bool const i_in_sensor = (mask[i] & groupbit);
    bool const i_in_source = !!(mask[i] & source_grpbit) != inv;
    if (!(i_in_sensor || i_in_source)) continue;
    double const xtmp = x[i][0];
    double const ytmp = x[i][1];
    double const ztmp = x[i][2];
    int const itype = type[i];
    bool const ipoint = !!ispoint[itype];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];
    for (int jj = 0; jj < jnum; jj++) {
      int const j = jlist[jj] & NEIGHMASK;
      bool const j_in_sensor = (mask[j] & groupbit);
      bool const j_in_source = !!(mask[j] & source_grpbit) != inv;
      bool const compute_ij = i_in_sensor && j_in_source;
      bool const compute_ji = (newton_pair || j < nlocal) && (j_in_sensor && i_in_source);
      if (!(compute_ij || compute_ji)) continue;
      double const delx = xtmp - x[j][0];
      double const dely = ytmp - x[j][1];
      double const delz = ztmp - x[j][2];
      double const rsq = delx * delx + dely * dely + delz * delz;
      int jtype = type[j];
      if (rsq >= cutsq[itype][jtype]) continue;
      double const factor_coul = special_coul[sbmask(j)];
      double const r = sqrt(rsq);
      double const rinv = 1.0 / r;
      double aij = rinv * ElectrodeMath::safe_erfc(alpha * r) - e_shift;
      double erfc_eta = 0.0;
      if (!(ipoint && !!ispoint[jtype])) {
        erfc_eta = ElectrodeMath::safe_erfc(eta[itype][jtype] * r);
        aij -= rinv * erfc_eta - eshift_eta[itype][jtype];
      }
      if (factor_coul < 1.0) aij -= (1.0 - factor_coul) * rinv * (1.0 - erfc_eta);
      if (i_in_sensor) { vec[i] += aij * q[j]; }
      if (j_in_sensor && (!inv || !i_in_sensor)) { vec[j] += aij * q[i]; }
    }
  }
}

/* ----------------------------------------------------------------------
   compute self interaction term of vector in constant potential method
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::compute_vector_self(double *vec, int groupbit, int source_grpbit,
                                                 bool inv)
{
  point_in_sensor_warning(groupbit);
  int const inum = list->inum;
  int *mask = atom->mask;
  int *type = atom->type;
  int *ilist = list->ilist;
  double *q = atom->q;

  double const selfint = 2.0 * alpha / MY_PIS;
  double const pre_eta = 2.0 / MY_PIS;
  double const pre_wolf = ElectrodeMath::safe_erfc(alpha * cut_coul) / cut_coul;

  for (int ii = 0; ii < inum; ii++) {
    int const i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    int const itype = type[i];
    bool const i_in_source = !!(mask[i] & source_grpbit) != inv;
    if (i_in_source) {
      vec[i] -= (selfint + pre_wolf) * q[i];
      if (!ispoint[itype])
        vec[i] += (pre_eta * eta[itype][itype] + eshift_eta[itype][itype]) * q[i];
    }
  }
}

/* ----------------------------------------------------------------------
   compute pair interaction term of matrix in constant potential method
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::compute_matrix(bigint *mpos, double **array, int groupbit)
{
  point_in_sensor_warning(groupbit);
  int *numneigh, **firstneigh;

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  int inum = list->inum;
  int *ilist = list->ilist;
  double *special_coul = force->special_coul;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double const e_shift = ElectrodeMath::safe_erfc(alpha * cut_coul) / cut_coul;

  for (int ii = 0; ii < inum; ii++) {
    int const i = ilist[ii];
    if (!(mask[i] & groupbit)) continue;
    int const itype = type[i];
    bool const ipoint = !!ispoint[itype];
    bigint const ipos = mpos[i];
    double const xtmp = x[i][0];
    double const ytmp = x[i][1];
    double const ztmp = x[i][2];
    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    // real-space part of matrix is symmetric
    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj];
      j &= NEIGHMASK;
      if (!(mask[j] & groupbit)) continue;

      double const delx = xtmp - x[j][0];    // neighlists take care of pbc
      double const dely = ytmp - x[j][1];
      double const delz = ztmp - x[j][2];
      double const rsq = delx * delx + dely * dely + delz * delz;
      int const jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        double const factor_coul = special_coul[sbmask(j)];
        double const r = sqrt(rsq);
        double const rinv = 1.0 / r;
        double aij = rinv * ElectrodeMath::safe_erfc(alpha * r) - e_shift;
        double erfc_eta = 0.0;
        if (!(ipoint && !!ispoint[jtype])) {
          double const erfc_eta = ElectrodeMath::safe_erfc(eta[itype][jtype] * r);
          aij -= rinv * erfc_eta - eshift_eta[itype][jtype];
        }
        if (factor_coul < 1.0) aij -= (1.0 - factor_coul) * rinv * (1.0 - erfc_eta);
        // newton on or off?
        if (!newton_pair && j >= nlocal) aij *= 0.5;
        bigint jpos = mpos[j];
        assert(jpos >= 0);
        array[ipos][jpos] += aij;
        array[jpos][ipos] += aij;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   compute self interaction term of matrix in constant potential method
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::compute_matrix_self(bigint *mpos, double **array, int groupbit)
{
  point_in_sensor_warning(groupbit);
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  int *type = atom->type;

  double const selfint = 2.0 * alpha / MY_PIS;
  double const pre_eta = 2.0 / MY_PIS;
  double const pre_wolf = ElectrodeMath::safe_erfc(alpha * cut_coul) / cut_coul;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      int const itype = type[i];
      array[mpos[i]][mpos[i]] -= selfint + pre_wolf;
      if (!ispoint[itype])
        array[mpos[i]][mpos[i]] += pre_eta * eta[itype][itype] + eshift_eta[itype][itype];
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        if (i == j) fwrite(&ispoint[i], sizeof(int), 1, fp);
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&eta[i][j], sizeof(double), 1, fp);
        fwrite(&cut_lj[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::read_restart(FILE *fp)
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
          if (i == j) utils::sfread(FLERR, &ispoint[i], sizeof(int), 1, fp, nullptr, error);
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &eta[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut_lj[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        if (i == j) MPI_Bcast(&ispoint[i], 1, MPI_INT, 0, world);
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&eta[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_lj[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::write_restart_settings(FILE *fp)
{
  fwrite(&alpha, sizeof(double), 1, fp);
  fwrite(&cut_lj_global, sizeof(double), 1, fp);
  fwrite(&cut_coul, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
  fwrite(&tabinner, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &alpha, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_lj_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_coul, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tabinner, sizeof(double), 1, fp, nullptr, error);
  }
  MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_lj_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_coul, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tabinner, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    if (!ispoint[i]) {
      fprintf(fp, "%d %g %g %g\n", i, epsilon[i][i], sigma[i][i], eta[i][i] * MY_SQRT2);
    } else {
      fprintf(fp, "%d %g %g %s\n", i, epsilon[i][i], sigma[i][i], "NULL");
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      if ((i == j) && !ispoint[i]) {
        fprintf(fp, "%d %d %g %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], eta[i][i] * MY_SQRT2,
                cut_lj[i][j]);
      } else {
        fprintf(fp, "%d %d %g %g %s %g\n", i, j, epsilon[i][j], sigma[i][j], "NULL", cut_lj[i][j]);
      }
}

/* ---------------------------------------------------------------------- */

void *PairLJCutCoulWolfGauss::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "alpha") == 0) return (void *) &alpha;
  if (strcmp(str, "cut_coul") == 0) return (void *) &cut_coul;
  dim = 1;
  if (strcmp(str, "ispoint") == 0) return (void *) ispoint;
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  if (strcmp(str, "eta") == 0) return (void *) eta;
  return nullptr;
}

/* ----------------------------------------------------------------------
   warn if point charges are in the sensor group (this should only be the case in EEM not for QEq or CPM)
------------------------------------------------------------------------- */

void PairLJCutCoulWolfGauss::point_in_sensor_warning(int groupbit)
{
  if (already_warned) return;
  int point_in_sensor = 0;
  for (int i = 0; i < atom->nlocal; i++) {
    if (!!ispoint[atom->type[i]] && (atom->mask[i] & groupbit)) point_in_sensor++;
  }
  MPI_Allreduce(MPI_IN_PLACE, &point_in_sensor, 1, MPI_INT, MPI_SUM, world);
  if (point_in_sensor) {
    if (comm->me == 0) error->warning(FLERR, "Point charges are used in sensor group");
    already_warned = true;
  }
}
