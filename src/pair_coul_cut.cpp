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

#include "pair_coul_cut.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulCut::PairCoulCut(LAMMPS *lmp) : Pair(lmp)
{
  born_matrix_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairCoulCut::~PairCoulCut()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulCut::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, ecoul, fpair;
  double rsq, r2inv, rinv, forcecoul, factor_coul;
  int *ilist, *jlist, *numneigh, **firstneigh;

  ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
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
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;
        rinv = sqrt(r2inv);
        forcecoul = qqrd2e * scale[itype][jtype] * qtmp * q[j] * rinv;
        fpair = factor_coul * forcecoul * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) ecoul = factor_coul * qqrd2e * scale[itype][jtype] * qtmp * q[j] * rinv;

        if (evflag) ev_tally(i, j, nlocal, newton_pair, 0.0, ecoul, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulCut::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(scale, np1, np1, "pair:scale");
  for (int i = 1; i < np1; i++) {
    for (int j = i; j < np1; j++) {
      setflag[i][j] = 0;
      scale[i][j] = 1.0;
    }
  }

  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(cut, np1, np1, "pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulCut::settings(int narg, char **arg)
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

void PairCoulCut::coeff(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double cut_one = cut_global;
  if (narg == 3) cut_one = utils::numeric(FLERR, arg[2], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      cut[i][j] = cut_one;
      scale[i][j] = 1.0;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulCut::init_style()
{
  if (!atom->q_flag) error->all(FLERR, "Pair style coul/cut requires atom attribute q");

  neighbor->add_request(this);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulCut::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i], cut[j][j]);
    scale[i][j] = 1.0;
  }

  scale[j][i] = scale[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCut::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&scale[i][j], sizeof(double), 1, fp);
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) { fwrite(&cut[i][j], sizeof(double), 1, fp); }
    }
  }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCut::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) {
        utils::sfread(FLERR, &scale[i][j], sizeof(double), 1, fp, nullptr, error);
        utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      }
      MPI_Bcast(&scale[i][j], 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulCut::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulCut::read_restart_settings(FILE *fp)
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

void PairCoulCut::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d\n", i);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairCoulCut::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) fprintf(fp, "%d %d %g\n", i, j, cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairCoulCut::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                           double factor_coul, double /*factor_lj*/, double &fforce)
{
  double r2inv, rinv, forcecoul, phicoul;

  r2inv = 1.0 / rsq;
  rinv = sqrt(r2inv);
  forcecoul = force->qqrd2e * atom->q[i] * atom->q[j] * rinv;
  fforce = factor_coul * forcecoul * r2inv;

  phicoul = force->qqrd2e * atom->q[i] * atom->q[j] * rinv;
  return factor_coul * phicoul;
}

/* ---------------------------------------------------------------------- */

void PairCoulCut::born_matrix(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                            double factor_coul, double /*factor_lj*/, double &dupair,
                            double &du2pair)
{
  double rinv, r2inv, r3inv;
  double du_coul, du2_coul;

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  r2inv = 1.0 / rsq;
  rinv = sqrt(r2inv);
  r3inv = r2inv * rinv;

  // Reminder: qqrd2e converts  q^2/r to energy w/ dielectric constant

  du_coul = -qqrd2e * q[i] * q[j] * r2inv;
  du2_coul = 2.0 * qqrd2e * q[i] * q[j] * r3inv;

  dupair = factor_coul * du_coul;
  du2pair = factor_coul * du2_coul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulCut::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "cut_coul") == 0) return (void *) cut;
  if (strcmp(str, "scale") == 0) return (void *) scale;
  return nullptr;
}
