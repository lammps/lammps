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

#include "pair_lj_expand_sphere.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_special.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathSpecial::powint;
using MathSpecial::square;

/* ---------------------------------------------------------------------- */

PairLJExpandSphere::PairLJExpandSphere(LAMMPS *lmp) :
    Pair(lmp), rmax(nullptr), cut(nullptr), epsilon(nullptr), sigma(nullptr), lj1(nullptr),
    lj2(nullptr), lj3(nullptr), lj4(nullptr)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJExpandSphere::~PairLJExpandSphere()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(rmax);
    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJExpandSphere::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, rtmp, delx, dely, delz, evdwl, fpair;
  double r, rshift, rshiftsq, rsq, r2inv, r6inv, forcelj, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *radius = atom->radius;
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
    rtmp = radius[i];
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

        // cutsq is maximum cutoff per type. Now compute and apply real cutoff

        r = sqrt(rsq);
        rshift = r - rtmp - radius[j];

        if (rshift < cut[itype][jtype]) {
          rshiftsq = rshift * rshift;

          r2inv = 1.0 / rshiftsq;
          r6inv = r2inv * r2inv * r2inv;

          forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
          fpair = factor_lj * forcelj * rshift / r;

          f[i][0] += delx * fpair;
          f[i][1] += dely * fpair;
          f[i][2] += delz * fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx * fpair;
            f[j][1] -= dely * fpair;
            f[j][2] -= delz * fpair;
          }

          if (eflag) {
            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]);
            if (offset_flag && (rshiftsq > 0.0)) {
              double ratio6 =
                  powint(sigma[itype][jtype] / (cut[itype][jtype] + rtmp + radius[j]), 6);
              evdwl -= 4.0 * epsilon[itype][jtype] * (ratio6 * ratio6 - ratio6);
            }
            evdwl *= factor_lj;
          }

          if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
        }
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJExpandSphere::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");

  memory->create(rmax, np1, "pair:rmax");
  memory->create(cut, np1, np1, "pair:cut");
  memory->create(epsilon, np1, np1, "pair:epsilon");
  memory->create(sigma, np1, np1, "pair:sigma");
  memory->create(lj1, np1, np1, "pair:lj1");
  memory->create(lj2, np1, np1, "pair:lj2");
  memory->create(lj3, np1, np1, "pair:lj3");
  memory->create(lj4, np1, np1, "pair:lj4");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJExpandSphere::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

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

void PairLJExpandSphere::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 5) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double cut_one = cut_global;
  if (narg == 5) cut_one = utils::numeric(FLERR, arg[4], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJExpandSphere::init_style()
{
  Pair::init_style();

  if (!atom->radius_flag)
    error->all(FLERR, "Pair style lj/expand/sphere requires atom attribute radius");

  // determine max radius per type

  int *type = atom->type;
  double *radius = atom->radius;
  rmax[0] = 0.0;
  for (int itype = 1; itype <= atom->ntypes; ++itype) {
    double rmax_one = 0.0;
    for (int i = 0; i < atom->nlocal; ++i) {
      if (type[i] == itype) rmax_one = MAX(rmax_one, radius[i]);
    }
    MPI_Allreduce(&rmax_one, &rmax[itype], 1, MPI_DOUBLE, MPI_MAX, world);
  }
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJExpandSphere::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], 0.0, 0.0);
    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
  }

  lj1[i][j] = 48.0 * epsilon[i][j] * powint(sigma[i][j], 12);
  lj2[i][j] = 24.0 * epsilon[i][j] * powint(sigma[i][j], 6);
  lj3[i][j] = 4.0 * epsilon[i][j] * powint(sigma[i][j], 12);
  lj4[i][j] = 4.0 * epsilon[i][j] * powint(sigma[i][j], 6);

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  cut[j][i] = cut[i][j];

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];

  // delta parameter is the sum of the radii
  return cut[i][j] + rmax[i] + rmax[j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJExpandSphere::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJExpandSphere::read_restart(FILE *fp)
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
          utils::sfread(FLERR, &epsilon[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &sigma[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJExpandSphere::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJExpandSphere::read_restart_settings(FILE *fp)
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

void PairLJExpandSphere::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d %g %g\n", i, epsilon[i][i], sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJExpandSphere::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJExpandSphere::single(int i, int j, int itype, int jtype, double rsq,
                                  double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r2inv, r6inv, r, rshift, rshiftsq, forcelj, philj;
  double *radius = atom->radius;

  fforce = philj = 0.0;
  r = sqrt(rsq);
  rshift = r - radius[i] - radius[j];
  if (rshift < cut[itype][jtype]) {
    rshiftsq = rshift * rshift;
    r2inv = 1.0 / rshiftsq;
    r6inv = r2inv * r2inv * r2inv;
    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
    fforce = factor_lj * forcelj * rshift / r;

    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]);
    if (offset_flag && (rshiftsq > 0.0)) {
      double ratio6 = powint(sigma[itype][jtype] / (cut[itype][jtype] + radius[i] + radius[j]), 6);
      philj -= 4.0 * epsilon[itype][jtype] * (ratio6 * ratio6 - ratio6);
    }
  }
  return factor_lj * philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJExpandSphere::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  return nullptr;
}
