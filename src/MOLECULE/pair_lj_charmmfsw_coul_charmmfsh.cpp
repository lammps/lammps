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
   Contributing author: Paul Crozier (SNL)
     The lj-fsw/coul-fsh (force-switched and force-shifted) sections
     were provided by Robert Meissner
     and Lucio Colombi Ciacchi of Bremen University, Bremen, Germany,
     with additional assistance from Robert A. Latour, Clemson University
------------------------------------------------------------------------- */

#include "pair_lj_charmmfsw_coul_charmmfsh.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "update.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCharmmfswCoulCharmmfsh::PairLJCharmmfswCoulCharmmfsh(LAMMPS *lmp) :
  Pair(lmp)
{
  implicit = 0;
  mix_flag = ARITHMETIC;
  writedata = 1;

  // short-range/long-range flag accessed by DihedralCharmmfsw

  dihedflag = 0;

  // switch qqr2e from LAMMPS value to CHARMM value

  if (strcmp(update->unit_style,"real") == 0) {
    if ((comm->me == 0) && (force->qqr2e != force->qqr2e_charmm_real))
      error->message(FLERR,"Switching to CHARMM coulomb energy"
                     " conversion constant");
    force->qqr2e = force->qqr2e_charmm_real;
  }
}

/* ---------------------------------------------------------------------- */

PairLJCharmmfswCoulCharmmfsh::~PairLJCharmmfswCoulCharmmfsh()
{
  // switch qqr2e back from CHARMM value to LAMMPS value

  if (update && strcmp(update->unit_style,"real") == 0) {
    if ((comm->me == 0) && (force->qqr2e == force->qqr2e_charmm_real))
      error->message(FLERR,"Restoring original LAMMPS coulomb energy"
                     " conversion constant");
    force->qqr2e = force->qqr2e_lammps_real;
  }

  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(eps14);
    memory->destroy(sigma14);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(lj14_1);
    memory->destroy(lj14_2);
    memory->destroy(lj14_3);
    memory->destroy(lj14_4);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,evdwl12,evdwl6,ecoul,fpair;
  double r,rinv,r3inv,rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double switch1;
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
        r = sqrt(rsq);

        if (rsq < cut_coulsq) {
          forcecoul = qqrd2e * qtmp*q[j]*
            (sqrt(r2inv) - r*cut_coulinv*cut_coulinv);
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          jtype = type[j];
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
              (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
            forcelj = forcelj*switch1;
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
            ecoul = qqrd2e * qtmp*q[j]*
              (sqrt(r2inv) + cut_coulinv*cut_coulinv*r - 2.0*cut_coulinv);
            ecoul *= factor_coul;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq) {
            if (rsq > cut_lj_innersq) {
              rinv = 1.0/r;
              r3inv = rinv*rinv*rinv;
              evdwl12 = lj3[itype][jtype]*cut_lj6*denom_lj12 *
                (r6inv - cut_lj6inv)*(r6inv - cut_lj6inv);
              evdwl6 = -lj4[itype][jtype]*cut_lj3*denom_lj6 *
                (r3inv - cut_lj3inv)*(r3inv - cut_lj3inv);;
              evdwl = evdwl12 + evdwl6;
            } else {
              evdwl12 = r6inv*lj3[itype][jtype]*r6inv -
                lj3[itype][jtype]*cut_lj_inner6inv*cut_lj6inv;
              evdwl6 = -lj4[itype][jtype]*r6inv +
                lj4[itype][jtype]*cut_lj_inner3inv*cut_lj3inv;
              evdwl = evdwl12 + evdwl6;
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

void PairLJCharmmfswCoulCharmmfsh::allocate()
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
  memory->create(eps14,n+1,n+1,"pair:eps14");
  memory->create(sigma14,n+1,n+1,"pair:sigma14");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(lj14_1,n+1,n+1,"pair:lj14_1");
  memory->create(lj14_2,n+1,n+1,"pair:lj14_2");
  memory->create(lj14_3,n+1,n+1,"pair:lj14_3");
  memory->create(lj14_4,n+1,n+1,"pair:lj14_4");
}

/* ----------------------------------------------------------------------
   global settings
   unlike other pair styles,
     there are no individual pair settings that these override
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::settings(int narg, char **arg)
{
  if (narg != 2 && narg != 3)
    error->all(FLERR,"Illegal pair_style command");

  cut_lj_inner = force->numeric(FLERR,arg[0]);
  cut_lj = force->numeric(FLERR,arg[1]);
  if (narg == 2) {
    cut_coul = cut_lj;
  } else {
    cut_coul = force->numeric(FLERR,arg[2]);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::coeff(int narg, char **arg)
{
  if (narg != 4 && narg != 6)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);
  double eps14_one = epsilon_one;
  double sigma14_one = sigma_one;
  if (narg == 6) {
    eps14_one = force->numeric(FLERR,arg[4]);
    sigma14_one = force->numeric(FLERR,arg[5]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      eps14[i][j] = eps14_one;
      sigma14[i][j] = sigma14_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/charmmfsw/coul/charmmfsh "
               "requires atom attribute q");

  neighbor->request(this,instance_me);

  // require cut_lj_inner < cut_lj

  if (cut_lj_inner >= cut_lj)
    error->all(FLERR,"Pair inner lj cutoff >= Pair outer lj cutoff");

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_ljinv = 1.0/cut_lj;
  cut_lj_innerinv = 1.0/cut_lj_inner;
  cut_lj3 = cut_lj * cut_lj * cut_lj;
  cut_lj3inv = cut_ljinv * cut_ljinv * cut_ljinv;
  cut_lj_inner3inv = cut_lj_innerinv * cut_lj_innerinv * cut_lj_innerinv;
  cut_lj_inner3 = cut_lj_inner * cut_lj_inner * cut_lj_inner;
  cut_lj6 = cut_ljsq * cut_ljsq * cut_ljsq;
  cut_lj6inv = cut_lj3inv * cut_lj3inv;
  cut_lj_inner6inv = cut_lj_inner3inv * cut_lj_inner3inv;
  cut_lj_inner6 = cut_lj_innersq * cut_lj_innersq * cut_lj_innersq;
  cut_coulsq = cut_coul * cut_coul;
  cut_coulinv = 1.0/cut_coul;
  cut_bothsq = MAX(cut_ljsq,cut_coulsq);

  denom_lj = (cut_ljsq-cut_lj_innersq) * (cut_ljsq-cut_lj_innersq) *
    (cut_ljsq-cut_lj_innersq);
  denom_lj12 = 1.0/(cut_lj6 - cut_lj_inner6);
  denom_lj6 = 1.0/(cut_lj3 - cut_lj_inner3);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCharmmfswCoulCharmmfsh::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    eps14[i][j] = mix_energy(eps14[i][i],eps14[j][j],
                               sigma14[i][i],sigma14[j][j]);
    sigma14[i][j] = mix_distance(sigma14[i][i],sigma14[j][j]);
  }

  double cut = MAX(cut_lj,cut_coul);

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj14_1[i][j] = 48.0 * eps14[i][j] * pow(sigma14[i][j],12.0);
  lj14_2[i][j] = 24.0 * eps14[i][j] * pow(sigma14[i][j],6.0);
  lj14_3[i][j] = 4.0 * eps14[i][j] * pow(sigma14[i][j],12.0);
  lj14_4[i][j] = 4.0 * eps14[i][j] * pow(sigma14[i][j],6.0);

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  lj14_1[j][i] = lj14_1[i][j];
  lj14_2[j][i] = lj14_2[i][j];
  lj14_3[j][i] = lj14_3[i][j];
  lj14_4[j][i] = lj14_4[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g\n",
            i,epsilon[i][i],sigma[i][i],eps14[i][i],sigma14[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,
              epsilon[i][j],sigma[i][j],eps14[i][j],sigma14[i][j]);
}


/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&eps14[i][j],sizeof(double),1,fp);
        fwrite(&sigma14[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&epsilon[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&eps14[i][j],sizeof(double),1,fp);
          fread(&sigma14[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&eps14[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma14[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::write_restart_settings(FILE *fp)
{
  fwrite(&cut_lj_inner,sizeof(double),1,fp);
  fwrite(&cut_lj,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCharmmfswCoulCharmmfsh::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_lj_inner,sizeof(double),1,fp);
    fread(&cut_lj,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_lj_inner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmfswCoulCharmmfsh::
single(int i, int j, int itype, int jtype,
       double rsq, double factor_coul, double factor_lj, double &fforce)
{
  double r,rinv,r2inv,r3inv,r6inv,forcecoul,forcelj;
  double phicoul,philj,philj12,philj6;
  double switch1;

  r2inv = 1.0/rsq;
  r = sqrt(rsq);
  rinv = 1.0/r;
  if (rsq < cut_coulsq) {
    forcecoul = force->qqrd2e * atom->q[i]*atom->q[j] *
      (sqrt(r2inv) - r*cut_coulinv*cut_coulinv);
  } else forcecoul = 0.0;

  if (rsq < cut_ljsq) {
    r6inv = r2inv*r2inv*r2inv;
    r3inv = rinv*rinv*rinv;
    forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
        (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      forcelj = forcelj*switch1;
    }
  } else forcelj = 0.0;

  fforce = (factor_coul*forcecoul + factor_lj*forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    phicoul = force->qqrd2e * atom->q[i]*atom->q[j] *
      (sqrt(r2inv) + cut_coulinv*cut_coulinv*r - 2.0*cut_coulinv);
    eng += factor_coul*phicoul;
  }
  if (rsq < cut_ljsq) {
    if (rsq > cut_lj_innersq) {
      philj12 = lj3[itype][jtype]*cut_lj6*denom_lj12 *
        (r6inv - cut_lj6inv)*(r6inv - cut_lj6inv);
      philj6 = -lj4[itype][jtype]*cut_lj3*denom_lj6 *
        (r3inv - cut_lj3inv)*(r3inv - cut_lj3inv);;
      philj = philj12 + philj6;
    } else {
      philj12 = r6inv*lj3[itype][jtype]*r6inv -
        lj3[itype][jtype]*cut_lj_inner6inv*cut_lj6inv;
      philj6 = -lj4[itype][jtype]*r6inv +
        lj4[itype][jtype]*cut_lj_inner3inv*cut_lj3inv;
      philj = philj12 + philj6;
    }
    eng += factor_lj*philj;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairLJCharmmfswCoulCharmmfsh::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"lj14_1") == 0) return (void *) lj14_1;
  if (strcmp(str,"lj14_2") == 0) return (void *) lj14_2;
  if (strcmp(str,"lj14_3") == 0) return (void *) lj14_3;
  if (strcmp(str,"lj14_4") == 0) return (void *) lj14_4;

  dim = 0;
  if (strcmp(str,"implicit") == 0) return (void *) &implicit;

  // info extracted by dihedral_charmmfsw

  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  if (strcmp(str,"cut_lj_inner") == 0) return (void *) &cut_lj_inner;
  if (strcmp(str,"cut_lj") == 0) return (void *) &cut_lj;
  if (strcmp(str,"dihedflag") == 0) return (void *) &dihedflag;

  return NULL;
}
