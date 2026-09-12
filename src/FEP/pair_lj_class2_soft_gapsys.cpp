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
   Adapted from the pair style lj/class2
------------------------------------------------------------------------- */

#include "pair_lj_class2_soft_gapsys.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJClass2SoftGapsys::PairLJClass2SoftGapsys(LAMMPS *lmp) :
    Pair(lmp), cut(nullptr), epsilon(nullptr), sigma(nullptr), lj1(nullptr), lj2(nullptr),
    lj3(nullptr), lj4(nullptr), lambda(nullptr), offset(nullptr)
{
  writedata = 1;
  allocated = 0;
  centroidstressflag = CENTROID_SAME;
}

/* ---------------------------------------------------------------------- */

PairLJClass2SoftGapsys::~PairLJClass2SoftGapsys()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lambda);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, rinv, r2inv, r3inv, r6inv, forcelj, factor_lj;
  double cut_inner;
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

      cut_inner = alphalj * pow(26.0 * lambda[itype][jtype] / 7.0, 1.0 / 6.0) * sigma[itype][jtype];

      if (rsq < cut_inner * cut_inner) {

        double epsln = epsilon[itype][jtype];

        double s_ri9 = pow(sigma[itype][jtype] / cut_inner, 9);
        double s_ri6 = pow(sigma[itype][jtype] / cut_inner, 6);

        double b1 = 18.0 * epsln * (7.0 * s_ri6 - 10.0 * s_ri9) / (cut_inner * cut_inner);
        double b2 = 18.0 * epsln * (11.0 * s_ri9 - 8.0 * s_ri6) / cut_inner;

        fpair = factor_lj * (b1 + b2 / sqrt(rsq));

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          double a3 = 110.0 * epsln * s_ri9 - 84 * epsln * s_ri6;
          evdwl = -0.5 * b1 * rsq - b2 * sqrt(rsq) + a3 - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);

      } else if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        rinv = sqrt(r2inv);
        r3inv = r2inv * rinv;
        r6inv = r3inv * r3inv;
        forcelj = r6inv * (lj1[itype][jtype] * r3inv - lj2[itype][jtype]);
        fpair = factor_lj * forcelj * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          evdwl = r6inv * (lj3[itype][jtype] * r3inv - lj4[itype][jtype]) - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::allocate()
{
  allocated = 1;
  const int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut, np1, np1, "pair:cut");
  memory->create(epsilon, np1, np1, "pair:epsilon");
  memory->create(sigma, np1, np1, "pair:sigma");
  memory->create(lambda, np1, np1, "pair:lambda");
  memory->create(lj1, np1, np1, "pair:lj1");
  memory->create(lj2, np1, np1, "pair:lj2");
  memory->create(lj3, np1, np1, "pair:lj3");
  memory->create(lj4, np1, np1, "pair:lj4");
  memory->create(offset, np1, np1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR, "Illegal pair_style command");

  alphalj = utils::numeric(FLERR, arg[0], false, lmp);
  if (!(alphalj > 0.0))
    error->all(FLERR, "Pair style lj/class2/soft/gapsys requires alphalj > 0");
  cut_global = utils::numeric(FLERR, arg[1], false, lmp);

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

void PairLJClass2SoftGapsys::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
  double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
  double lambda_one = utils::numeric(FLERR,arg[4],false,lmp);

  if (sigma_one <= 0.0 || epsilon_one <= 0.0)
    error->all(FLERR,"Incorrect args for pair coefficients" + utils::errorurl(21));
  if (lambda_one < 0.0 || lambda_one > 1.0)
    error->all(FLERR, "Pair style lj/class2/soft/gapsys requires 0 <= lambda <= 1");

  double cut_one = cut_global;
  if (narg == 6) cut_one = utils::numeric(FLERR, arg[5], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      lambda[i][j] = lambda_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients" + utils::errorurl(21));
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJClass2SoftGapsys::init_one(int i, int j)
{
  // always mix epsilon,sigma via sixthpower rules
  // mix distance via user-defined rule

  if (setflag[i][j] == 0) {
    epsilon[i][j] = 2.0 * sqrt(epsilon[i][i] * epsilon[j][j]) * pow(sigma[i][i], 3.0) *
        pow(sigma[j][j], 3.0) / (pow(sigma[i][i], 6.0) + pow(sigma[j][j], 6.0));
    sigma[i][j] = pow((0.5 * (pow(sigma[i][i], 6.0) + pow(sigma[j][j], 6.0))), 1.0 / 6.0);
    if (lambda[i][i] != lambda[j][j])
      error->all(FLERR,"Pair lj/class2/soft/gapsys different lambda values in mix");
    lambda[i][j] = lambda[i][i];
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
    did_mix = true;
  }

  lj1[i][j] = 18.0 * epsilon[i][j] * pow(sigma[i][j], 9.0);
  lj2[i][j] = 18.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  lj3[i][j] = 2.0 * epsilon[i][j] * pow(sigma[i][j], 9.0);
  lj4[i][j] = 3.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

  if (offset_flag && (cut[i][j] > 0.0)) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = epsilon[i][j] * (2.0 * pow(ratio, 9.0) - 3.0 * pow(ratio, 6.0));
  } else
    offset[i][j] = 0.0;

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  lambda[j][i] = lambda[i][j];
  cut[j][i] = cut[i][j];
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

    double sig3 = sigma[i][j] * sigma[i][j] * sigma[i][j];
    double sig6 = sig3 * sig3;
    double rc3 = cut[i][j] * cut[i][j] * cut[i][j];
    double rc6 = rc3 * rc3;
    etail_ij =
        2.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 * (sig3 - 3.0 * rc3) / (3.0 * rc6);
    ptail_ij = 2.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 * (sig3 - 2.0 * rc3) / rc6;
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j], sizeof(double), 1, fp);
        fwrite(&sigma[i][j], sizeof(double), 1, fp);
        fwrite(&lambda[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::read_restart(FILE *fp)
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
          utils::sfread(FLERR, &lambda[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&lambda[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::write_restart_settings(FILE *fp)
{
  fwrite(&alphalj, sizeof(double), 1, fp);
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
  fwrite(&tail_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&alphalj, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tail_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&alphalj, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&tail_flag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, epsilon[i][i], sigma[i][i], lambda[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJClass2SoftGapsys::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g\n", i, j, epsilon[i][j], sigma[i][j], lambda[i][j], cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJClass2SoftGapsys::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                            double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r2inv, rinv, r3inv, r6inv, forcelj, philj;

  double cut_inner;
  cut_inner = alphalj * pow(26.0 * lambda[itype][jtype] / 7.0, 1.0 / 6.0) * sigma[itype][jtype];

  if (rsq > cut_inner * cut_inner) {

    r2inv = 1.0 / rsq;
    rinv = sqrt(r2inv);
    r3inv = r2inv * rinv;
    r6inv = r3inv * r3inv;
    forcelj = r6inv * (lj1[itype][jtype] * r3inv - lj2[itype][jtype]);
    fforce = factor_lj * forcelj * r2inv;

    philj = r6inv * (lj3[itype][jtype] * r3inv - lj4[itype][jtype]) - offset[itype][jtype];
    return factor_lj * philj;

  } else {

    double epsln = epsilon[itype][jtype];

    double s_ri9 = pow(sigma[itype][jtype] / cut_inner, 9);
    double s_ri6 = pow(sigma[itype][jtype] / cut_inner, 6);

    double b1 = 18.0 * epsln * (7.0 * s_ri6 - 10.0 * s_ri9) / (cut_inner * cut_inner);
    double b2 = 18.0 * epsln * (11.0 * s_ri9 - 8.0 * s_ri6) / cut_inner;
    fforce = factor_lj * (b1 + b2 / sqrt(rsq));

    double a3 = 110.0 * epsln * s_ri9 - 84 * epsln * s_ri6;
    return factor_lj * (-0.5 * b1 * rsq - b2 * sqrt(rsq) + a3 - offset[itype][jtype]);

  }
}

/* ---------------------------------------------------------------------- */

void *PairLJClass2SoftGapsys::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "epsilon") == 0) return (void *) epsilon;
  if (strcmp(str, "sigma") == 0) return (void *) sigma;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  return nullptr;
}
