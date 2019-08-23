/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ----------------------------------------------------------------------
   Contributing authors: Eugen Rozic (University College London)
------------------------------------------------------------------------- */

#include "pair_cosine_squared.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairCosineSquared::PairCosineSquared(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairCosineSquared::~PairCosineSquared()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(w);
    memory->destroy(cut);
    memory->destroy(wcaflag);

    memory->destroy(lj12_e);
    memory->destroy(lj6_e);
    memory->destroy(lj12_f);
    memory->destroy(lj6_f);
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCosineSquared::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n+1, n+1, "pair:setflag");
  memory->create(cutsq, n+1, n+1, "pair:cutsq");

  memory->create(cut, n+1, n+1, "pair:cut");
  memory->create(epsilon, n+1, n+1, "pair:epsilon");
  memory->create(sigma, n+1, n+1, "pair:sigma");
  memory->create(w, n+1, n+1, "pair:w");
  memory->create(wcaflag, n+1, n+1, "pair:wcaflag");

  memory->create(lj12_e, n+1, n+1, "pair:lj12_e");
  memory->create(lj6_e, n+1, n+1, "pair:lj6_e");
  memory->create(lj12_f, n+1, n+1, "pair:lj12_f");
  memory->create(lj6_f, n+1, n+1, "pair:lj6_f");

  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      setflag[i][j] = 0;
      wcaflag[i][j] = 0;
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCosineSquared::settings(int narg, char **arg)
{
  if (narg != 1) {
    error->all(FLERR, "Illegal pair_style command (wrong number of params)");
  }

  cut_global = force->numeric(FLERR, arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j])
          cut[i][j] = cut_global;
  }
}


/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCosineSquared::coeff(int narg, char **arg)
{
  if (narg < 4 || narg > 6)
    error->all(FLERR, "Incorrect args for pair coefficients (too few or too many)");
  
  if (!allocated)
    allocate();

  int ilo, ihi, jlo, jhi;
  force->bounds(FLERR, arg[0], atom->ntypes, ilo, ihi);
  force->bounds(FLERR, arg[1], atom->ntypes, jlo, jhi);

  double epsilon_one = force->numeric(FLERR, arg[2]);
  double sigma_one = force->numeric(FLERR, arg[3]);

  double cut_one = cut_global;
  double wca_one = 0;
  if (narg == 6) {
    cut_one = force->numeric(FLERR, arg[4]);
    if (strcmp(arg[5], "wca") == 0) {
      wca_one = 1;
    } else {
      error->all(FLERR, "Incorrect args for pair coefficients (unknown option)");
    }
  } else if (narg == 5) {
    if (strcmp(arg[4], "wca") == 0) {
      wca_one = 1;
    } else {
      cut_one = force->numeric(FLERR, arg[4]);
    }
  }

  if (cut_one < sigma_one) {
    error->all(FLERR, "Incorrect args for pair coefficients (cutoff < sigma)");
  } else if (cut_one == sigma_one) {
    if (wca_one == 0) {
      error->all(FLERR, "Incorrect args for pair coefficients (cutoff = sigma w/o wca)");
    } else {
      error->warning(FLERR, "Cosine/squared set to WCA only (cutoff = sigma)");
    }
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      cut[i][j] = cut_one;
      wcaflag[i][j] = wca_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0)
    error->all(FLERR, "Incorrect args for pair coefficients (none set)");
}

/* ----------------------------------------------------------------------
   init specific to this pair style (unneccesary)
------------------------------------------------------------------------- */

/*
void PairCosineSquared::init_style()
{
  neighbor->request(this,instance_me);
}
*/

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCosineSquared::init_one(int i, int j)
{
  if (setflag[i][j] == 0)
    error->all(FLERR, "Mixing not supported in pair_style cosine/squared");

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  cut[j][i] = cut[i][j];
  wcaflag[j][i] = wcaflag[i][j];

  w[j][i] = w[i][j] = cut[i][j] - sigma[i][j];

  if (wcaflag[i][j]) {
    lj12_e[j][i] = lj12_e[i][j] = epsilon[i][j] * pow(sigma[i][j], 12.0);
    lj6_e[j][i] = lj6_e[i][j] = 2.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
    lj12_f[j][i] = lj12_f[i][j] = 12.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
    lj6_f[j][i] = lj6_f[i][j] = 12.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
  }

  // Note: cutsq is set in pair.cpp

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   this is here to throw errors & warnings for given options
------------------------------------------------------------------------- */

void PairCosineSquared::modify_params(int narg, char **arg)
{
  Pair::modify_params(narg, arg);

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "mix") == 0) {
      error->all(FLERR, "pair_modify mix not supported for pair_style cosine/squared");
    } else if (strcmp(arg[iarg], "shift") == 0) {
      error->warning(FLERR, "pair_modify shift has no effect on pair_style cosine/squared");
      offset_flag = 0;
    } else if (strcmp(arg[iarg], "tail") == 0) {
      error->warning(FLERR, "pair_modify tail has no effect on pair_style cosine/squared");
      tail_flag = 0;
    }
    iarg++;
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCosineSquared::write_restart(FILE *fp)
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
        fwrite(&wcaflag[i][j], sizeof(int), 1, fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCosineSquared::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0)
        fread(&setflag[i][j], sizeof(int), 1, fp);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j], sizeof(double), 1, fp);
          fread(&sigma[i][j], sizeof(double), 1, fp);
          fread(&cut[i][j], sizeof(double), 1, fp);
          fread(&wcaflag[i][j], sizeof(int), 1, fp);
        }
        MPI_Bcast(&epsilon[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&sigma[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&wcaflag[i][j], 1, MPI_INT, 0, world);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCosineSquared::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global, sizeof(double), 1, fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCosineSquared::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global, sizeof(double), 1, fp);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairCosineSquared::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g %d\n", i, epsilon[i][i], sigma[i][i],
            cut[i][i], wcaflag[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairCosineSquared::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %d\n", i, j, epsilon[i][j], sigma[i][j],
              cut[i][j], wcaflag[i][j]);
}

/* ---------------------------------------------------------------------- */

void PairCosineSquared::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double r, rsq, r2inv, r6inv;
  double factor_lj, force_lj, force_cos, cosone;

  evdwl = 0.0;
  if (eflag || vflag)
    ev_setup(eflag, vflag);
  else
    evflag = vflag_fdotr = 0;

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
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {

        /*
          This is exactly what the "single" method does, in fact it could be called
          here instead of repeating the code but here energy calculation is optional
          so a little bit of calculation is possibly saved
        */

        r = sqrt(rsq);

        if (r <= sigma[itype][jtype]) {
          if (wcaflag[itype][jtype]) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            force_lj = r6inv*(lj12_f[itype][jtype]*r6inv - lj6_f[itype][jtype]);
            fpair = factor_lj*force_lj*r2inv;
            if (eflag) {
              evdwl = factor_lj*r6inv*
                      (lj12_e[itype][jtype]*r6inv - lj6_e[itype][jtype]);
              if (sigma[itype][jtype] == cut[itype][jtype]) {
                // this is the WCA-only case (it requires this shift by definition)
                evdwl += factor_lj*epsilon[itype][jtype];
              }
            }
          } else {
            fpair = 0.0;
            if (eflag) {
              evdwl = -factor_lj*epsilon[itype][jtype];
            }
          }
        } else {
          force_cos = -(MY_PI*epsilon[itype][jtype] / (2.0*w[itype][jtype])) *
                      sin(MY_PI*(r-sigma[itype][jtype]) / w[itype][jtype]);
          fpair = factor_lj*force_cos / r;
          if (eflag) {
            cosone = cos(MY_PI*(r-sigma[itype][jtype]) / (2.0*w[itype][jtype]));
            evdwl = -factor_lj*epsilon[itype][jtype]*cosone*cosone;
          }
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;

        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr)
    virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   This is used be pair_write;
   it is called only if rsq < cutsq[itype][jtype], no need to check that
------------------------------------------------------------------------- */

double PairCosineSquared::single(int /* i */, int /* j */, int itype, int jtype, double rsq,
                         double /* factor_coul */, double factor_lj,
                         double &fforce)
{
  double r, r2inv, r6inv, cosone, force, energy;
  
  r = sqrt(rsq);

  if (r <= sigma[itype][jtype]) {
    if (wcaflag[itype][jtype]) {
      r2inv = 1.0/rsq;
      r6inv = r2inv*r2inv*r2inv;
      force = r6inv*(lj12_f[itype][jtype]*r6inv - lj6_f[itype][jtype])*r2inv;
      energy = r6inv*(lj12_e[itype][jtype]*r6inv - lj6_e[itype][jtype]);
      if (sigma[itype][jtype] == cut[itype][jtype]) {
        // this is the WCA-only case (it requires this shift by definition)
        energy += epsilon[itype][jtype];
      }
    } else {
      force = 0.0;
      energy = -epsilon[itype][jtype];
    }
  } else {
    cosone = cos(MY_PI*(r-sigma[itype][jtype]) / (2.0*w[itype][jtype]));
    force = -(MY_PI*epsilon[itype][jtype] / (2.0*w[itype][jtype])) * 
                 sin(MY_PI*(r-sigma[itype][jtype]) / w[itype][jtype]) / r;
    energy = -epsilon[itype][jtype]*cosone*cosone;
  }
  fforce = factor_lj*force;
  return factor_lj*energy;
}

