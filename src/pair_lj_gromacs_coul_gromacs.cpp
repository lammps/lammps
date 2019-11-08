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
   Contributing author: Mark Stevens (SNL)
------------------------------------------------------------------------- */

#include "pair_lj_gromacs_coul_gromacs.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJGromacsCoulGromacs::PairLJGromacsCoulGromacs(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJGromacsCoulGromacs::~PairLJGromacsCoulGromacs()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(ljsw1);
    memory->destroy(ljsw2);
    memory->destroy(ljsw3);
    memory->destroy(ljsw4);
    memory->destroy(ljsw5);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double r,tlj,tc,fswitch,fswitchcoul,eswitch,ecoulswitch;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

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
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_bothsq) {
        r2inv = 1.0/rsq;

        // skip if qi or qj = 0.0 since this potential may be used as
        // coarse-grain model with many uncharged atoms

        if (rsq < cut_coulsq && qtmp != 0.0 && q[j] != 0.0) {
          forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
          if (rsq > cut_coul_innersq) {
            r = sqrt(rsq);
            tc = r - cut_coul_inner;
            fswitchcoul = qqrd2e * qtmp*q[j]*r*tc*tc*(coulsw1 + coulsw2*tc);
            forcecoul += fswitchcoul;
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          jtype = type[j];
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          if (rsq > cut_lj_innersq) {
            r = sqrt(rsq);
            tlj = r - cut_lj_inner;
            fswitch = r*tlj*tlj*(ljsw1[itype][jtype] +
                                 ljsw2[itype][jtype]*tlj);
            forcelj += fswitch;
          }
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

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = qqrd2e * qtmp*q[j] * (sqrt(r2inv) - coulsw5);
            if (rsq > cut_coul_innersq) {
              ecoulswitch = tc*tc*tc * (coulsw3 + coulsw4*tc);
              ecoul += qqrd2e*qtmp*q[j]*ecoulswitch;
            }
            ecoul *= factor_coul;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq) {
            evdwl = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
            evdwl += ljsw5[itype][jtype];
            if (rsq > cut_lj_innersq) {
              eswitch = tlj*tlj*tlj *
                (ljsw3[itype][jtype] + ljsw4[itype][jtype]*tlj);
              evdwl += eswitch;
            }
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(ljsw1,n+1,n+1,"pair:ljsw1");
  memory->create(ljsw2,n+1,n+1,"pair:ljsw2");
  memory->create(ljsw3,n+1,n+1,"pair:ljsw3");
  memory->create(ljsw4,n+1,n+1,"pair:ljsw4");
  memory->create(ljsw5,n+1,n+1,"pair:ljsw5");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::settings(int narg, char **arg)
{
  if (narg != 2 && narg != 4)
    error->all(FLERR,"Illegal pair_style command");

  cut_lj_inner = force->numeric(FLERR,arg[0]);
  cut_lj = force->numeric(FLERR,arg[1]);
  if (narg == 2) {
    cut_coul_inner = cut_lj_inner;
    cut_coul = cut_lj;
  } else {
    cut_coul_inner = force->numeric(FLERR,arg[2]);
    cut_coul = force->numeric(FLERR,arg[3]);
  }

  if (cut_lj_inner <= 0.0 || cut_coul_inner < 0.0)
    error->all(FLERR,"Illegal pair_style command");
  if (cut_lj_inner > cut_lj || cut_coul_inner > cut_coul)
    error->all(FLERR,"Illegal pair_style command");
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,
               "Pair style lj/gromacs/coul/gromacs requires atom attribute q");

  neighbor->request(this,instance_me);

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coul_innersq = cut_coul_inner * cut_coul_inner;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq,cut_coulsq);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJGromacsCoulGromacs::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
  }

  double cut = MAX(cut_lj,cut_coul);

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);

  double r6inv = 1.0/pow(cut_lj,6.0);
  double r8inv = 1.0/pow(cut_lj,8.0);
  double t = cut_lj - cut_lj_inner;
  double t2inv = 1.0/(t*t);
  double t3inv = t2inv/t;
  double t3 = 1.0/t3inv;
  double a6 = (7.0*cut_lj_inner - 10.0*cut_lj)*r8inv*t2inv;
  double b6 = (9.0*cut_lj -  7.0*cut_lj_inner)*r8inv*t3inv;
  double a12 = (13.0*cut_lj_inner - 16.0*cut_lj)*r6inv*r8inv*t2inv;
  double b12 = (15.0*cut_lj - 13.0*cut_lj_inner)*r6inv*r8inv*t3inv;
  double c6 = r6inv - t3*(6.0*a6/3.0 + 6.0*b6*t/4.0);
  double c12 = r6inv*r6inv - t3*(12.0*a12/3.0 + 12.0*b12*t/4.0);

  ljsw1[i][j] = lj1[i][j]*a12 - lj2[i][j]*a6;
  ljsw2[i][j] = lj1[i][j]*b12 - lj2[i][j]*b6;
  ljsw3[i][j] = -lj3[i][j]*12.0*a12/3.0 + lj4[i][j]*6.0*a6/3.0;
  ljsw4[i][j] = -lj3[i][j]*12.0*b12/4.0 + lj4[i][j]*6.0*b6/4.0;
  ljsw5[i][j] = -lj3[i][j]*c12 + lj4[i][j]*c6;

  double r3inv = 1.0/pow(cut_coul,3.0);
  t = cut_coul - cut_coul_inner;
  t2inv = 1.0/(t*t);
  t3inv = t2inv/t;
  double a1 = (2.0*cut_coul_inner - 5.0*cut_coul) * r3inv*t2inv;
  double b1 = (4.0*cut_coul - 2.0*cut_coul_inner) * r3inv*t3inv;
  coulsw1 = a1;
  coulsw2 = b1;
  coulsw3 = -a1/3.0;
  coulsw4 = -b1/4.0;
  coulsw5 = 1.0/cut_coul - t*t*t*(a1/3.0 + b1*t/4.0);

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  ljsw1[j][i] = ljsw1[i][j];
  ljsw2[j][i] = ljsw2[i][j];
  ljsw3[j][i] = ljsw3[i][j];
  ljsw4[j][i] = ljsw4[i][j];
  ljsw5[j][i] = ljsw5[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,NULL,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&epsilon[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&sigma[i][j],sizeof(double),1,fp,NULL,error);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_inner,sizeof(double),1,fp);
  fwrite(&cut_lj,sizeof(double),1,fp);
  fwrite(&cut_coul_inner,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&cut_lj_inner,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_lj,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_coul_inner,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&cut_lj_inner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul_inner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJGromacsCoulGromacs::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g\n",i,j,epsilon[i][j],sigma[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJGromacsCoulGromacs::single(int i, int j, int itype, int jtype,
                                double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,forcecoul,forcelj,phicoul,philj;
  double r,tlj,tc,fswitch,phiswitch,fswitchcoul,phiswitchcoul;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j]*sqrt(r2inv);
    if (rsq > cut_coul_innersq) {
      r = sqrt(rsq);
      tc = r - cut_coul_inner;
      fswitchcoul =  force->qqrd2e *
        atom->q[i]*atom->q[j] * r*tc*tc * (coulsw1 + coulsw2*tc);
      forcecoul += fswitchcoul;
    }
  } else forcecoul = 0.0;

  if (rsq < cut_ljsq) {
    r6inv = r2inv*r2inv*r2inv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    if (rsq > cut_lj_innersq) {
      r = sqrt(rsq);
      tlj = r - cut_lj_inner;
      fswitch = r*tlj*tlj*(ljsw1[itype][jtype] + ljsw2[itype][jtype]*tlj);
      forcelj += fswitch;
    }
  } else forcelj = 0.0;

  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j] * (sqrt(r2inv)-coulsw5);
    if (rsq > cut_coul_innersq) {
      phiswitchcoul = force->qqrd2e * atom->q[i]*atom->q[j] *
        tc*tc*tc * (coulsw3 + coulsw4*tc);
      phicoul += phiswitchcoul;
    }
    eng += factor_coul*phicoul;
  }

  if (rsq < cut_ljsq) {
    philj = r6inv * (lj3[itype][jtype]*r6inv - lj4[itype][jtype]);
    philj += ljsw5[itype][jtype];
    if (rsq > cut_lj_innersq) {
      phiswitch = tlj*tlj*tlj *
        (ljsw3[itype][jtype] + ljsw4[itype][jtype]*tlj);
      philj += phiswitch;
    }
    eng += factor_lj*philj;
  }

  return eng;
}
