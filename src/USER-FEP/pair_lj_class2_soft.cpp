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

#include "pair_lj_class2_soft.h"
#include <mpi.h>
#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJClass2Soft::PairLJClass2Soft(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;
  cntratmstressflag = 1;
}

/* ---------------------------------------------------------------------- */

PairLJClass2Soft::~PairLJClass2Soft()
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
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJClass2Soft::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,forcelj,factor_lj;
  double denlj, r4sig6;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

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
        denlj = lj3[itype][jtype] + pow(rsq, 3) * pow(sigma[itype][jtype], -6.0);
        r4sig6 = rsq*rsq / lj2[itype][jtype];
        forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
            (18.0*r4sig6*pow(denlj, -2.5) - 18.0*r4sig6*pow(denlj, -2));
        fpair = factor_lj*forcelj;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          denlj = lj3[itype][jtype] + pow(rsq, 3) * pow(sigma[itype][jtype], -6.0);
          evdwl = lj1[itype][jtype] * epsilon[itype][jtype] * (2.0/(denlj*sqrt(denlj)) - 3.0/denlj) -
            offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJClass2Soft::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(epsilon,n+1,n+1,"pair:epsilon");
  memory->create(sigma,n+1,n+1,"pair:sigma");
  memory->create(lambda,n+1,n+1,"pair:lambda");
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJClass2Soft::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  nlambda = force->numeric(FLERR,arg[0]);
  alphalj = force->numeric(FLERR,arg[1]);

  cut_global = force->numeric(FLERR,arg[2]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJClass2Soft::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);
  double lambda_one = force->numeric(FLERR,arg[4]);
  if (sigma_one <= 0.0) error->all(FLERR,"Incorrect args for pair coefficients");

  double cut_one = cut_global;
  if (narg == 6) cut_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      lambda[i][j] = lambda_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJClass2Soft::init_one(int i, int j)
{
  // always mix epsilon,sigma via sixthpower rules
  // mix distance via user-defined rule

  if (setflag[i][j] == 0) {
    epsilon[i][j] = 2.0 * sqrt(epsilon[i][i]*epsilon[j][j]) *
      pow(sigma[i][i],3.0) * pow(sigma[j][j],3.0) /
      (pow(sigma[i][i],6.0) + pow(sigma[j][j],6.0));
    sigma[i][j] =
      pow((0.5 * (pow(sigma[i][i],6.0) + pow(sigma[j][j],6.0))),1.0/6.0);
    if (lambda[i][i] != lambda[j][j])
      error->all(FLERR,"Pair lj/class2/coul/cut/soft different lambda values in mix");
    lambda[i][j] = lambda[i][i];
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  lj1[i][j] = pow(lambda[i][j], nlambda);
  lj2[i][j] = pow(sigma[i][j], 6.0);
  lj3[i][j] = alphalj * (1.0 - lambda[i][j])*(1.0 - lambda[i][j]);

  if (offset_flag && (cut[i][j] > 0.0)) {
    double denlj = lj3[i][j] + pow(cut[i][j] / sigma[i][j], 6.0);
    offset[i][j] = lj1[i][j] * epsilon[i][j] * (2.0/(denlj*sqrt(denlj)) - 3.0/denlj);
  } else offset[i][j] = 0.0;

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  lambda[j][i] = lambda[i][j];
  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  offset[j][i] = offset[i][j];

  // compute I,J contribution to long-range tail correction
  // count total # of atoms of type I and J via Allreduce

  if (tail_flag) {
    int *type = atom->type;
    int nlocal = atom->nlocal;

    double count[2],all[2];
    count[0] = count[1] = 0.0;
    for (int k = 0; k < nlocal; k++) {
      if (type[k] == i) count[0] += 1.0;
      if (type[k] == j) count[1] += 1.0;
    }
    MPI_Allreduce(count,all,2,MPI_DOUBLE,MPI_SUM,world);

    double sig3 = sigma[i][j]*sigma[i][j]*sigma[i][j];
    double sig6 = sig3*sig3;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    // check the following expressions for etail_lj & ptail_ij they are not correct !
    etail_ij = 2.0*MY_PI*all[0]*all[1]*lj1[i][j] *epsilon[i][j] *
      sig6 * (sig3 - 3.0*rc3) / (3.0*rc6);
    ptail_ij = 2.0*MY_PI*all[0]*all[1]*lj1[i][j] *epsilon[i][j] *
      sig6 * (sig3 - 2.0*rc3) / rc6;
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJClass2Soft::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&lambda[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJClass2Soft::read_restart(FILE *fp)
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
          utils::sfread(FLERR,&lambda[i][j],sizeof(double),1,fp,NULL,error);
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,NULL,error);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lambda[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJClass2Soft::write_restart_settings(FILE *fp)
{
  fwrite(&nlambda,sizeof(double),1,fp);
  fwrite(&alphalj,sizeof(double),1,fp);
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJClass2Soft::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&nlambda,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&alphalj,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,NULL,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,NULL,error);
    utils::sfread(FLERR,&tail_flag,sizeof(int),1,fp,NULL,error);
  }
  MPI_Bcast(&nlambda,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alphalj,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJClass2Soft::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g\n",i,epsilon[i][i],sigma[i][i],lambda[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJClass2Soft::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],
              lambda[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJClass2Soft::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                            double /*factor_coul*/, double factor_lj,
                            double &fforce)
{
  double forcelj,philj;
  double r4sig6, denlj;

  if (rsq < cutsq[itype][jtype]) {
    r4sig6 = rsq*rsq / lj2[itype][jtype];
    denlj = lj3[itype][jtype] + rsq*r4sig6;
    forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
      (18.0*r4sig6/(denlj*denlj*sqrt(denlj)) - 18.0*r4sig6/(denlj*denlj));
  } else forcelj = 0.0;
  fforce = factor_lj*forcelj;

  if (rsq < cutsq[itype][jtype]) {
    denlj = lj3[itype][jtype] + rsq*rsq*rsq / lj2[itype][jtype];
    philj = lj1[itype][jtype] * epsilon[itype][jtype] * (2.0/(denlj*sqrt(denlj)) - 3.0/denlj) -
      offset[itype][jtype];
  } else philj = 0.0;

  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJClass2Soft::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;
  return NULL;
}
