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
   Contributing author: Carsten Svaneborg (SDU)
------------------------------------------------------------------------- */

#include "pair_zero2.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "memory.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairZero2::PairZero2(LAMMPS *lmp) : Pair(lmp)
{
  coeffflag = 1;
  writedata = 1;
  single_enable = 1;
  respa_enable = 1;
}

/* ---------------------------------------------------------------------- */

PairZero2::~PairZero2()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
  }
}

/* ---------------------------------------------------------------------- */

void PairZero2::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);
  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairZero2::compute_outer(int eflag, int vflag)
{
  ev_init(eflag, vflag);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairZero2::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) setflag[i][j] = 0;

  memory->create(cutsq, n + 1, n + 1, "pair:cutsq");
  memory->create(cut, n + 1, n + 1, "pair:cut");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairZero2::settings(int narg, char **arg)
{
  if ((narg != 1) && (narg != 2)) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (narg == 2) {
    if (strcmp("nocoeff", arg[1]) == 0)
      coeffflag = 0;
    else
      error->all(FLERR, "Illegal pair_style command");
  }

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i + 1; j <= atom->ntypes; j++) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairZero2::coeff(int narg, char **arg)
{
  if ((narg < 2) || (coeffflag && narg > 3))
    error->all(FLERR, "Incorrect args for pair coefficients");

  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double cut_one = cut_global;
  if (coeffflag && (narg == 3)) cut_one = utils::numeric(FLERR, arg[2], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
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

double PairZero2::init_one(int i, int j)
{
  if (setflag[i][j] == 0) { cut[i][j] = mix_distance(cut[i][i], cut[j][j]); }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairZero2::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) { fwrite(&cut[i][j], sizeof(double), 1, fp); }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairZero2::read_restart(FILE *fp)
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
        if (me == 0) { utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error); }
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairZero2::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&coeffflag, sizeof(int), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairZero2::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &coeffflag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&coeffflag, 1, MPI_INT, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairZero2::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++) fprintf(fp, "%d\n", i);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairZero2::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++) fprintf(fp, "%d %d %g\n", i, j, cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairZero2::single(int /*i*/, int /*j*/, int /* itype */, int /* jtype */, double /* rsq */,
                         double /*factor_coul*/, double /* factor_lj */, double &fforce)
{
  fforce = 0.0;
  return 0.0;
}
