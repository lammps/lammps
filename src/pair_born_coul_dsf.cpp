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
   Contributing author: Ariel Lozano (arielzn@gmail.com)
   References: Fennell and Gezelter, JCP 124, 234104 (2006)
------------------------------------------------------------------------- */

#include "pair_born_coul_dsf.h"
#include <mpi.h>
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairBornCoulDSF::PairBornCoulDSF(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  single_enable = 0;
}

/* ---------------------------------------------------------------------- */

PairBornCoulDSF::~PairBornCoulDSF()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut_lj);
    memory->destroy(cut_ljsq);
    memory->destroy(a);
    memory->destroy(rho);
    memory->destroy(sigma);
    memory->destroy(c);
    memory->destroy(d);
    memory->destroy(rhoinv);
    memory->destroy(born1);
    memory->destroy(born2);
    memory->destroy(born3);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairBornCoulDSF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,r2inv,r6inv,forcecoul,forceborn,factor_coul,factor_lj;
  double prefactor,erfcc,erfcd,arg;
  double rexp;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
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

    // self coulombic energy

    if (eflag) {
      double e_self = -(e_shift/2.0 + alpha/MY_PIS) * qtmp*qtmp*qqrd2e;
      ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          prefactor = qqrd2e*qtmp*q[j]/r;
          arg = alpha * r ;
          erfcd = MathSpecial::expmsq(arg);
          erfcc = MathSpecial::my_erfcx(arg) * erfcd;
          forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd +
                                   r*f_shift) * r;
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          r = sqrt(rsq);
          rexp = exp((sigma[itype][jtype]-r)*rhoinv[itype][jtype]);
          forceborn = born1[itype][jtype]*r*rexp - born2[itype][jtype]*r6inv
            + born3[itype][jtype]*r2inv*r6inv;
        } else forceborn = 0.0;

        fpair = (forcecoul + factor_lj*forceborn) * r2inv;

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
            ecoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv +
              d[itype][jtype]*r6inv*r2inv - offset[itype][jtype];
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

void PairBornCoulDSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut_lj,n+1,n+1,"pair:cut_lj");
  memory->create(cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  memory->create(a,n+1,n+1,"pair:a");
  memory->create(rho,n+1,n+1,"pair:rho");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(c,n+1,n+1,"pair:c");
  memory->create(d,n+1,n+1,"pair:d");
  memory->create(rhoinv,n+1,n+1,"pair:rhoinv");
  memory->create(born1,n+1,n+1,"pair:born1");
  memory->create(born2,n+1,n+1,"pair:born2");
  memory->create(born3,n+1,n+1,"pair:born3");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBornCoulDSF::settings(int narg, char **arg)
{
  if (narg < 2 || narg > 3) error->all(FLERR,"Illegal pair_style command");

  alpha = force->numeric(FLERR,arg[0]);
  cut_lj_global = force->numeric(FLERR,arg[1]);
  if (narg == 2) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(FLERR,arg[2]);

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBornCoulDSF::coeff(int narg, char **arg)
{
  if (narg < 7 || narg > 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double a_one = force->numeric(FLERR,arg[2]);
  double rho_one = force->numeric(FLERR,arg[3]);
  double sigma_one = force->numeric(FLERR,arg[4]);
  if (rho_one <= 0) error->all(FLERR,"Incorrect args for pair coefficients");
  double c_one = force->numeric(FLERR,arg[5]);
  double d_one = force->numeric(FLERR,arg[6]);

  double cut_lj_one = cut_lj_global;
  if (narg == 8) cut_lj_one = force->numeric(FLERR,arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      a[i][j] = a_one;
      rho[i][j] = rho_one;
      sigma[i][j] = sigma_one;
      c[i][j] = c_one;
      d[i][j] = d_one;
      cut_lj[i][j] = cut_lj_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairBornCoulDSF::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style born/coul/dsf requires atom attribute q");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;
  double erfcc = erfc(alpha*cut_coul);
  double erfcd = exp(-alpha*alpha*cut_coul*cut_coul);
  f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*alpha*erfcd/cut_coul);
  e_shift = erfcc/cut_coul - f_shift*cut_coul;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairBornCoulDSF::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  double cut = MAX(cut_lj[i][j],cut_coul);
  cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];

  rhoinv[i][j] = 1.0/rho[i][j];
  born1[i][j] = a[i][j]/rho[i][j];
  born2[i][j] = 6.0*c[i][j];
  born3[i][j] = 8.0*d[i][j];

  if (offset_flag && (cut_lj[i][j] > 0.0)) {
    double rexp = exp((sigma[i][j]-cut_lj[i][j])*rhoinv[i][j]);
    offset[i][j] = a[i][j]*rexp - c[i][j]/pow(cut_lj[i][j],6.0)
      + d[i][j]/pow(cut_lj[i][j],8.0);
  } else offset[i][j] = 0.0;

  cut_ljsq[j][i] = cut_ljsq[i][j];
  a[j][i] = a[i][j];
  c[j][i] = c[i][j];
  d[j][i] = d[i][j];
  rhoinv[j][i] = rhoinv[i][j];
  sigma[j][i] = sigma[i][j];
  born1[j][i] = born1[i][j];
  born2[j][i] = born2[i][j];
  born3[j][i] = born3[i][j];
  offset[j][i] = offset[i][j];

  return cut;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBornCoulDSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&a[i][j],sizeof(double),1,fp);
        fwrite(&rho[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&c[i][j],sizeof(double),1,fp);
        fwrite(&d[i][j],sizeof(double),1,fp);
        fwrite(&cut_lj[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBornCoulDSF::read_restart(FILE *fp)
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
          fread(&a[i][j],sizeof(double),1,fp);
          fread(&rho[i][j],sizeof(double),1,fp);
          fread(&sigma[i][j],sizeof(double),1,fp);
          fread(&c[i][j],sizeof(double),1,fp);
          fread(&d[i][j],sizeof(double),1,fp);
          fread(&cut_lj[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&a[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&rho[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&c[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&d[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut_lj[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBornCoulDSF::write_restart_settings(FILE *fp)
{
  fwrite(&alpha,sizeof(double),1,fp);
  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBornCoulDSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&alpha,sizeof(double),1,fp);
    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&alpha,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairBornCoulDSF::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,
            a[i][i],rho[i][i],sigma[i][i],c[i][i],d[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairBornCoulDSF::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g %g\n",i,j,
              a[i][j],rho[i][j],sigma[i][j],c[i][j],d[i][j],cut_lj[i][j]);
}

/* ----------------------------------------------------------------------
   only the pair part is calculated here
------------------------------------------------------------------------- */

double PairBornCoulDSF::single(int i, int j, int itype, int jtype,
                                double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,r,prefactor,rexp;
  double forcecoul,forceborn,phicoul,phiborn;
  double erfcc,erfcd,arg;

  r2inv = 1.0/rsq;

  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = factor_coul * force->qqrd2e * atom->q[i]*atom->q[j]/r;

    arg = alpha * r ;
    erfcd = MathSpecial::expmsq(arg);
    erfcc = MathSpecial::my_erfcx(arg) * erfcd;

    forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS*erfcd +
      r*f_shift) * r;

    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  } else forcecoul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    r = sqrt(rsq);
    rexp = exp(-r*rhoinv[itype][jtype]);
    forceborn = born1[itype][jtype]*r*rexp - born2[itype][jtype]*r6inv
      + born3[itype][jtype]*r2inv*r6inv;
  } else forceborn = 0.0;

  fforce = (forcecoul + factor_lj*forceborn) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    phicoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    eng += phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    phiborn = a[itype][jtype]*rexp - c[itype][jtype]*r6inv
      + d[itype][jtype]*r2inv*r6inv - offset[itype][jtype];
    eng += factor_lj*phiborn;
  }
  return eng;
}
