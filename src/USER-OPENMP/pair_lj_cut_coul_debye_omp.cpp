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

#include "math.h"
#include "stdlib.h"
#include "pair_lj_cut_coul_debye_omp.h"
#include "atom.h"
#include "neigh_list.h"
#include "force.h"
#include "comm.h"
#include "error.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutCoulDebyeOMP::PairLJCutCoulDebyeOMP(LAMMPS *lmp) : PairLJCutCoulCutOMP(lmp) {}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulDebyeOMP::compute(int eflag, int vflag)
{
  if (eflag || vflag) {
    ev_setup(eflag,vflag);
    ev_setup_thr(eflag,vflag);
  } else evflag = vflag_fdotr = 0;

  if (evflag) {
    if (eflag) {
      if (force->newton_pair) return eval<1,1,1>();
      else return eval<1,1,0>();
    } else {
      if (force->newton_pair) return eval<1,0,1>();
      else return eval<1,0,0>();
    }
  } else {
    if (force->newton_pair) return eval<0,0,1>();
    else return eval<0,0,0>();
  }
}

template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
void PairLJCutCoulDebyeOMP::eval()
{

#if defined(_OPENMP)
#pragma omp parallel default(shared)
#endif
    {
      int i,j,ii,jj,inum,jnum,itype,jtype,tid;
      double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
      double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
      double r,rinv,screening;
      int *ilist,*jlist,*numneigh,**firstneigh;

      evdwl = ecoul = 0.0;

      const int nlocal = atom->nlocal;
      const int nall = nlocal + atom->nghost;
      const int nthreads = comm->nthreads;

      double **x = atom->x;
      double **f = atom->f;
      double *q = atom->q;
      int *type = atom->type;
      double *special_coul = force->special_coul;
      double *special_lj = force->special_lj;
      int newton_pair = force->newton_pair;
      double qqrd2e = force->qqrd2e;

      inum = list->inum;
      ilist = list->ilist;
      numneigh = list->numneigh;
      firstneigh = list->firstneigh;

      // loop over neighbors of my atoms

      int iifrom, iito;
      f = loop_setup_thr(f, iifrom, iito, tid, inum, nall, nthreads);
      for (ii = iifrom; ii < iito; ++ii) {
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

              if (j < nall) factor_coul = factor_lj = 1.0;
              else {
                  factor_coul = special_coul[j/nall];
                  factor_lj = special_lj[j/nall];
                  j %= nall;
              }

              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx*delx + dely*dely + delz*delz;
              jtype = type[j];

              if (rsq < cutsq[itype][jtype]) {
                  r2inv = 1.0/rsq;

                  if (rsq < cut_coulsq[itype][jtype]) {
                      r = sqrt(rsq);
                      rinv = 1.0/r;
                      screening = exp(-kappa*r);
                      forcecoul = qqrd2e * qtmp*q[j] * screening * (kappa + rinv);
                  } else forcecoul = 0.0;

                  if (rsq < cut_ljsq[itype][jtype]) {
                      r6inv = r2inv*r2inv*r2inv;
                      forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
                  } else forcelj = 0.0;

                  fpair = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

                  f[i][0] += delx*fpair;
                  f[i][1] += dely*fpair;
                  f[i][2] += delz*fpair;
                  if (newton_pair || j < nlocal) {
                      f[j][0] -= delx*fpair;
                      f[j][1] -= dely*fpair;
                      f[j][2] -= delz*fpair;
                  }

                  if (EFLAG) {
                  if (rsq < cut_coulsq[itype][jtype])
                    ecoul = factor_coul * qqrd2e * qtmp*q[j] * rinv * screening;
                  else ecoul = 0.0;
                      if (rsq < cut_ljsq[itype][jtype]) {
                        evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
                          offset[itype][jtype];
                        evdwl *= factor_lj;
                      } else evdwl = 0.0;
                  }

                  if (EVFLAG) ev_tally(i,j,nlocal,newton_pair,
                      evdwl,ecoul,fpair,delx,dely,delz);
              }
          }
      }
      // reduce per thread forces into global force array.
      force_reduce_thr(atom->f, nall, nthreads, tid);
    }
    ev_reduce_thr();
    if (vflag_fdotr) virial_compute();
}

/* ----------------------------------------------------------------------
   global settings 
------------------------------------------------------------------------- */

void PairLJCutCoulDebyeOMP::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all("Illegal pair_style command");

  kappa = force->numeric(arg[0]);
  cut_lj_global = force->numeric(arg[1]);
  if (narg == 2) cut_coul_global = cut_lj_global;
  else cut_coul_global = force->numeric(arg[2]);

  // reset cutoffs that were previously set from data file

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j] == 1) {
	  cut_lj[i][j] = cut_lj_global;
	  cut_coul[i][j] = cut_coul_global;
	}
  }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulDebyeOMP::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul_global,sizeof(double),1,fp);
  fwrite(&kappa,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulDebyeOMP::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul_global,sizeof(double),1,fp);
    fread(&kappa,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&kappa,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulDebyeOMP::single(int i, int j, int itype, int jtype,
				  double rsq,
				  double factor_coul, double factor_lj,
				  double &fforce)
{
  double r2inv,r6inv,r,rinv,screening,forcecoul,forcelj,phicoul,philj;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq[itype][jtype]) {
    r = sqrt(rsq);
    rinv = 1.0/r;
    screening = exp(-kappa*r);
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] *
      screening * (kappa + rinv);
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  } else forcelj = 0.0;
  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j] * rinv * screening;
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
      offset[itype][jtype];
    eng += factor_lj*philj;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulDebyeOMP::memory_usage()
{
  const int n=atom->ntypes;

  double bytes = PairOMP::memory_usage();

  bytes += 9*((n+1)*(n+1) * sizeof(double) + (n+1)*sizeof(double *));
  bytes += 1*((n+1)*(n+1) * sizeof(int) + (n+1)*sizeof(int *));

  return bytes;
}
