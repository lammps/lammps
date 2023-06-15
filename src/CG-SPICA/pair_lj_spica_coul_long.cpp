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
   Contributing author: Axel Kohlmeyer (Temple U)
   This style is a simplified re-implementation of the CG/CMM pair style
------------------------------------------------------------------------- */

#include "pair_lj_spica_coul_long.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

#define LMP_NEED_SPICA_FIND_LJ_TYPE 1
#include "lj_spica_common.h"

using namespace LAMMPS_NS;
using namespace LJSPICAParms;

#define EWALD_F 1.12837917
#define EWALD_P 0.3275911
#define A1 0.254829592
#define A2 -0.284496736
#define A3 1.421413741
#define A4 -1.453152027
#define A5 1.061405429

/* ---------------------------------------------------------------------- */

PairLJSPICACoulLong::PairLJSPICACoulLong(LAMMPS *lmp) :
    Pair(lmp), cut_lj(nullptr), cut_ljsq(nullptr), epsilon(nullptr), sigma(nullptr), lj1(nullptr),
    lj2(nullptr), lj3(nullptr), lj4(nullptr), offset(nullptr), lj_type(nullptr), rminsq(nullptr),
    emin(nullptr)
{
  ewaldflag = pppmflag = 1;
  respa_enable = 0;
  writedata = 1;
  ftable = nullptr;
}

/* ---------------------------------------------------------------------- */

PairLJSPICACoulLong::~PairLJSPICACoulLong()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(lj_type);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);

    memory->destroy(rminsq);
    memory->destroy(emin);

    allocated = 0;
  }
  if (ftable) free_tables();
}

/* ---------------------------------------------------------------------- */

void PairLJSPICACoulLong::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (evflag) {
    if (eflag) {
      if (force->newton_pair)
        eval<1, 1, 1>();
      else
        eval<1, 1, 0>();
    } else {
      if (force->newton_pair)
        eval<1, 0, 1>();
      else
        eval<1, 0, 0>();
    }
  } else {
    if (force->newton_pair)
      eval<0, 0, 1>();
    else
      eval<0, 0, 0>();
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void PairLJSPICACoulLong::eval()
{
  int i, ii, j, jj, jtype, itable;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double fraction, table;
  double r, rsq, r2inv, forcecoul, forcelj, factor_coul, factor_lj;
  double grij, expm2, prefactor, t, erfc;

  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const double *const q = atom->q;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  const double *const special_coul = force->special_coul;
  const double *const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;
  double fxtmp, fytmp, fztmp;

  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {

    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    const int itype = type[i];
    const int *const jlist = firstneigh[i];
    const int jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      forcecoul = forcelj = evdwl = ecoul = 0.0;

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
        const int ljt = lj_type[itype][jtype];

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij * grij);
            t = 1.0 / (1.0 + EWALD_P * grij);
            erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
            prefactor = qqrd2e * qtmp * q[j] / r;
            forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
            if (EFLAG) ecoul = prefactor * erfc;
            if (factor_coul < 1.0) {
              forcecoul -= (1.0 - factor_coul) * prefactor;
              if (EFLAG) ecoul -= (1.0 - factor_coul) * prefactor;
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction * dftable[itable];
            forcecoul = qtmp * q[j] * table;
            if (EFLAG) ecoul = qtmp * q[j] * (etable[itable] + fraction * detable[itable]);
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction * dctable[itable];
              prefactor = qtmp * q[j] * table;
              forcecoul -= (1.0 - factor_coul) * prefactor;
              if (EFLAG) ecoul -= (1.0 - factor_coul) * prefactor;
            }
          }
        }

        if (rsq < cut_ljsq[itype][jtype]) {

          if (ljt == LJ12_4) {
            const double r4inv = r2inv * r2inv;
            forcelj = r4inv * (lj1[itype][jtype] * r4inv * r4inv - lj2[itype][jtype]);

            if (EFLAG)
              evdwl = r4inv * (lj3[itype][jtype] * r4inv * r4inv - lj4[itype][jtype]) -
                  offset[itype][jtype];

          } else if (ljt == LJ9_6) {
            const double r3inv = r2inv * sqrt(r2inv);
            const double r6inv = r3inv * r3inv;
            forcelj = r6inv * (lj1[itype][jtype] * r3inv - lj2[itype][jtype]);
            if (EFLAG)
              evdwl =
                  r6inv * (lj3[itype][jtype] * r3inv - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_6) {
            const double r6inv = r2inv * r2inv * r2inv;
            forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
            if (EFLAG)
              evdwl =
                  r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];

          } else if (ljt == LJ12_5) {
            const double r5inv = r2inv * r2inv * sqrt(r2inv);
            const double r7inv = r5inv * r2inv;
            forcelj = r5inv * (lj1[itype][jtype] * r7inv - lj2[itype][jtype]);
            if (EFLAG)
              evdwl =
                  r5inv * (lj3[itype][jtype] * r7inv - lj4[itype][jtype]) - offset[itype][jtype];
          }
          forcelj *= factor_lj;
          if (EFLAG) evdwl *= factor_lj;
        }

        fpair = (forcecoul + forcelj) * r2inv;

        fxtmp += delx * fpair;
        fytmp += dely * fpair;
        fztmp += delz * fpair;
        if (NEWTON_PAIR || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (EVFLAG) ev_tally(i, j, nlocal, NEWTON_PAIR, evdwl, ecoul, fpair, delx, dely, delz);
      }
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(lj_type, np1, np1, "pair:lj_type");
  for (int i = 1; i < np1; i++) {
    for (int j = i; j < np1; j++) {
      setflag[i][j] = 0;
      lj_type[i][j] = LJ_NOT_SET;
    }
  }

  memory->create(cutsq, np1, np1, "pair:cutsq");

  memory->create(cut_lj, np1, np1, "pair:cut_lj");
  memory->create(cut_ljsq, np1, np1, "pair:cut_ljsq");
  memory->create(epsilon, np1, np1, "pair:epsilon");
  memory->create(sigma, np1, np1, "pair:sigma");

  memory->create(lj1, np1, np1, "pair:lj1");
  memory->create(lj2, np1, np1, "pair:lj2");
  memory->create(lj3, np1, np1, "pair:lj3");
  memory->create(lj4, np1, np1, "pair:lj4");

  memory->create(offset, np1, np1, "pair:offset");

  memory->create(rminsq, np1, np1, "pair:rminsq");
  memory->create(emin, np1, np1, "pair:emin");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR, "Illegal pair_style command");

  cut_lj_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (narg == 1)
    cut_coul = cut_lj_global;
  else
    cut_coul = utils::numeric(FLERR, arg[1], false, lmp);

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

void PairLJSPICACoulLong::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  int lj_type_one = find_lj_type(arg[2], lj_type_list);
  if (lj_type_one == LJ_NOT_SET) error->all(FLERR, "Cannot parse LJ type flag.");

  double epsilon_one = utils::numeric(FLERR, arg[3], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[4], false, lmp);

  double cut_lj_one = cut_lj_global;
  if (narg == 6) cut_lj_one = utils::numeric(FLERR, arg[5], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      lj_type[i][j] = lj_type_one;
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style lj/cut/coul/long requires atom attribute q");

  neighbor->add_request(this);

  cut_coulsq = cut_coul * cut_coul;

  // ensure use of KSpace long-range solver, set g_ewald

  if (force->kspace == nullptr) error->all(FLERR, "Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;

  // setup force tables (no rRESPA support yet)

  if (ncoultablebits) init_tables(cut_coul, nullptr);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJSPICACoulLong::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR,
               "No mixing support for lj/spica/coul/long. "
               "Coefficients for all pairs need to be set explicitly.");

  const int ljt = lj_type[i][j];

  if (ljt == LJ_NOT_SET) error->all(FLERR, "unrecognized LJ parameter flag");

  double cut = MAX(cut_lj[i][j], cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  lj1[i][j] = lj_prefact[ljt] * lj_pow1[ljt] * epsilon[i][j] * pow(sigma[i][j], lj_pow1[ljt]);
  lj2[i][j] = lj_prefact[ljt] * lj_pow2[ljt] * epsilon[i][j] * pow(sigma[i][j], lj_pow2[ljt]);
  lj3[i][j] = lj_prefact[ljt] * epsilon[i][j] * pow(sigma[i][j], lj_pow1[ljt]);
  lj4[i][j] = lj_prefact[ljt] * epsilon[i][j] * pow(sigma[i][j], lj_pow2[ljt]);

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut_lj[i][j];
    offset[i][j] =
        lj_prefact[ljt] * epsilon[i][j] * (pow(ratio, lj_pow1[ljt]) - pow(ratio, lj_pow2[ljt]));
  } else
    offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  cut_lj[j][i] = cut_lj[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];
  lj_type[j][i] = lj_type[i][j];

  // compute LJ derived parameters for SPICA angle potential (LJ only!)

  const double eps = epsilon[i][j];
  const double sig = sigma[i][j];
  const double rmin =
      sig * exp(1.0 / (lj_pow1[ljt] - lj_pow2[ljt]) * log(lj_pow1[ljt] / lj_pow2[ljt]));
  rminsq[j][i] = rminsq[i][j] = rmin * rmin;

  const double ratio = sig / rmin;
  const double emin_one =
      lj_prefact[ljt] * eps * (pow(ratio, lj_pow1[ljt]) - pow(ratio, lj_pow2[ljt]));
  emin[j][i] = emin[i][j] = emin_one;

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) error->all(FLERR, "Tail flag not supported by lj/spica/coul/long pair style");

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&lj_type[i][j], sizeof(int), 1, fp);
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut_lj[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::read_restart(FILE *fp)
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
          utils::sfread(FLERR, &lj_type[i][j], sizeof(int), 1, fp, nullptr, error);
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut_lj[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&lj_type[i][j], 1, MPI_INT, 0, world);
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut_lj[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global, sizeof(double), 1, fp);
  fwrite(&cut_coul, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
  fwrite(&ncoultablebits, sizeof(int), 1, fp);
  fwrite(&tabinner, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_lj_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_coul, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &ncoultablebits, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tabinner, sizeof(double), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_lj_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_coul, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&ncoultablebits, 1, MPI_INT, 0, world);
  MPI_Bcast(&tabinner, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   lj/spica does not support per atom type output with mixing
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::write_data(FILE *)
{
  error->one(FLERR,
             "Pair style lj/spica/coul/* requires using "
             "write_data with the 'pair ij' option");
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJSPICACoulLong::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %s %g %g %g\n", i, j, lj_type_list[lj_type[i][j]], epsilon[i][j],
              sigma[i][j], cut_lj[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJSPICACoulLong::single(int i, int j, int itype, int jtype, double rsq,
                                   double factor_coul, double factor_lj, double &fforce)
{
  double r2inv, r, grij, expm2, t, erfc, prefactor;
  double fraction, table, forcecoul, forcelj, phicoul, philj;
  int itable;

  forcecoul = forcelj = phicoul = philj = 0.0;

  r2inv = 1.0 / rsq;
  if (rsq < cut_coulsq) {
    if (!ncoultablebits || rsq <= tabinnersq) {
      r = sqrt(rsq);
      grij = g_ewald * r;
      expm2 = exp(-grij * grij);
      t = 1.0 / (1.0 + EWALD_P * grij);
      erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;
      prefactor = force->qqrd2e * atom->q[i] * atom->q[j] / r;
      forcecoul = prefactor * (erfc + EWALD_F * grij * expm2);
      phicoul = prefactor * erfc;
      if (factor_coul < 1.0) {
        forcecoul -= (1.0 - factor_coul) * prefactor;
        phicoul -= (1.0 - factor_coul) * prefactor;
      }
    } else {
      union_int_float_t rsq_lookup_single;
      rsq_lookup_single.f = rsq;
      itable = rsq_lookup_single.i & ncoulmask;
      itable >>= ncoulshiftbits;
      fraction = (rsq_lookup_single.f - rtable[itable]) * drtable[itable];
      table = ftable[itable] + fraction * dftable[itable];
      forcecoul = atom->q[i] * atom->q[j] * table;
      table = etable[itable] + fraction * detable[itable];
      phicoul = atom->q[i] * atom->q[j] * table;
      if (factor_coul < 1.0) {
        table = ctable[itable] + fraction * dctable[itable];
        prefactor = atom->q[i] * atom->q[j] * table;
        forcecoul -= (1.0 - factor_coul) * prefactor;
        phicoul -= (1.0 - factor_coul) * prefactor;
      }
    }
  }

  if (rsq < cut_ljsq[itype][jtype]) {
    const int ljt = lj_type[itype][jtype];

    if (ljt == LJ12_4) {
      const double r4inv = r2inv * r2inv;
      forcelj = r4inv * (lj1[itype][jtype] * r4inv * r4inv - lj2[itype][jtype]);

      philj =
          r4inv * (lj3[itype][jtype] * r4inv * r4inv - lj4[itype][jtype]) - offset[itype][jtype];

    } else if (ljt == LJ9_6) {
      const double r3inv = r2inv * sqrt(r2inv);
      const double r6inv = r3inv * r3inv;
      forcelj = r6inv * (lj1[itype][jtype] * r3inv - lj2[itype][jtype]);
      philj = r6inv * (lj3[itype][jtype] * r3inv - lj4[itype][jtype]) - offset[itype][jtype];

    } else if (ljt == LJ12_6) {
      const double r6inv = r2inv * r2inv * r2inv;
      forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
      philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];

    } else if (ljt == LJ12_5) {
      const double r5inv = r2inv * r2inv * sqrt(r2inv);
      const double r7inv = r5inv * r2inv;
      forcelj = r5inv * (lj1[itype][jtype] * r7inv - lj2[itype][jtype]);
      philj = r5inv * (lj3[itype][jtype] * r7inv - lj4[itype][jtype]) - offset[itype][jtype];
    }
    forcelj *= factor_lj;
    philj *= factor_lj;
  }

  fforce = (forcecoul + forcelj) * r2inv;

  return phicoul + philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJSPICACoulLong::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  if (strcmp(str, "lj_type") == 0) return (void *) lj_type;
  if (strcmp(str, "lj1") == 0) return (void *) lj1;
  if (strcmp(str, "lj2") == 0) return (void *) lj2;
  if (strcmp(str, "lj3") == 0) return (void *) lj3;
  if (strcmp(str, "lj4") == 0) return (void *) lj4;
  if (strcmp(str, "rminsq") == 0) return (void *) rminsq;
  if (strcmp(str, "emin") == 0) return (void *) emin;

  dim = 0;
  if (strcmp(str, "cut_coul") == 0) return (void *) &cut_coul;
  return nullptr;
}

/* ---------------------------------------------------------------------- */

double PairLJSPICACoulLong::memory_usage()
{
  double bytes = Pair::memory_usage();
  int n = atom->ntypes;

  // setflag/lj_type
  bytes += (double) 2 * (n + 1) * (n + 1) * sizeof(int);
  // lj_cut/lj_cutsq/epsilon/sigma/offset/lj1/lj2/lj3/lj4/rminsq/emin
  bytes += (double) 11 * (n + 1) * (n + 1) * sizeof(double);

  if (ncoultablebits) {
    int ntable = 1 << ncoultablebits;
    bytes += (double) 8 * ntable * sizeof(double);
  }

  return bytes;
}
