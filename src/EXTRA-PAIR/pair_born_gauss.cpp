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

// Contributing author: Axel Kohlmeyer, Temple University, akohlmey@gmail.com

#include "pair_born_gauss.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairBornGauss::PairBornGauss(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 1;
  respa_enable = 0;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairBornGauss::~PairBornGauss()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(biga0);
    memory->destroy(alpha);
    memory->destroy(biga1);
    memory->destroy(beta);
    memory->destroy(r0);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairBornGauss::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r, dr, aexp, bexp, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

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
        dr = r - r0[itype][jtype];
        aexp = biga0[itype][jtype] * exp(-alpha[itype][jtype] * r);
        bexp = biga1[itype][jtype] * exp(-beta[itype][jtype] * dr * dr);
        fpair = alpha[itype][jtype] * aexp;
        fpair -= 2.0 * beta[itype][jtype] * dr * bexp;
        fpair *= factor_lj / r;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) evdwl = factor_lj * (aexp - bexp - offset[itype][jtype]);
        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairBornGauss::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut, np1, np1, "pair:cut");
  memory->create(biga0, np1, np1, "pair:biga0");
  memory->create(alpha, np1, np1, "pair:alpha");
  memory->create(biga1, np1, np1, "pair:biga1");
  memory->create(beta, np1, np1, "pair:beta");
  memory->create(r0, np1, np1, "pair:r0");
  memory->create(offset, np1, np1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBornGauss::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Pair style bond/gauss must have exactly one argument");
  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset per-type pair cutoffs that have been explicitly set previously

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBornGauss::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double biga0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double alpha_one = utils::numeric(FLERR, arg[3], false, lmp);
  double biga1_one = utils::numeric(FLERR, arg[4], false, lmp);
  double beta_one = utils::numeric(FLERR, arg[5], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[6], false, lmp);
  double cut_one = cut_global;
  if (narg == 10) cut_one = utils::numeric(FLERR, arg[7], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      biga0[i][j] = biga0_one;
      alpha[i][j] = alpha_one;
      biga1[i][j] = biga1_one;
      beta[i][j] = beta_one;
      r0[i][j] = r0_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBornGauss::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  if (offset_flag) {
    double dr = cut[i][j] - r0[i][j];
    offset[i][j] =
        biga0[i][j] * exp(-alpha[i][j] * cut[i][j]) - biga1[i][j] * exp(-beta[i][j] * dr * dr);
  } else
    offset[i][j] = 0.0;

  biga0[j][i] = biga0[i][j];
  alpha[j][i] = alpha[i][j];
  biga1[j][i] = biga1[i][j];
  beta[j][i] = beta[i][j];
  r0[j][i] = r0[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBornGauss::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&biga0[i][j], sizeof(double), 1, fp);
        fwrite(&alpha[i][j], sizeof(double), 1, fp);
        fwrite(&biga1[i][j], sizeof(double), 1, fp);
        fwrite(&beta[i][j], sizeof(double), 1, fp);
        fwrite(&r0[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBornGauss::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &biga0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &alpha[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &biga1[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &beta[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &r0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&biga0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&alpha[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&biga1[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&beta[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&r0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBornGauss::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBornGauss::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
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

void PairBornGauss::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %g %g\n", i, biga0[i][i], alpha[i][i], biga1[i][i], beta[i][i],
            r0[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBornGauss::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g %g %g\n", i, j, biga0[i][j], alpha[i][j], biga1[i][j],
              beta[i][j], r0[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairBornGauss::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                             double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r, dr, aexp, bexp;

  r = sqrt(rsq);
  dr = r - r0[itype][jtype];
  aexp = biga0[itype][jtype] * exp(-alpha[itype][jtype] * r);
  bexp = biga1[itype][jtype] * exp(-beta[itype][jtype] * dr * dr);

  fforce = factor_lj * (alpha[itype][jtype] * aexp - 2.0 * dr * beta[itype][jtype] * bexp) / r;
  return factor_lj * (aexp - bexp - offset[itype][jtype]);
}

/* ---------------------------------------------------------------------- */

void *PairBornGauss::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "biga0") == 0) return (void *) biga0;
  if (strcmp(str, "biga1") == 0) return (void *) biga1;
  if (strcmp(str, "r0") == 0) return (void *) r0;
  return nullptr;
}
