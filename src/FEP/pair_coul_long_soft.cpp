// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
   Soft-core version: Agilio Padua (ENS de Lyon & CNRS)
------------------------------------------------------------------------- */

#include "pair_coul_long_soft.h"

#include <cmath>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"


using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairCoulLongSoft::PairCoulLongSoft(LAMMPS *lmp) : Pair(lmp)
{
  ewaldflag = pppmflag = 1;
  qdist = 0.0;
}

/* ---------------------------------------------------------------------- */

PairCoulLongSoft::~PairCoulLongSoft()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(scale);

    memory->destroy(lambda);
    memory->destroy(lam1);
    memory->destroy(lam2);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulLongSoft::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,rsq,forcecoul,factor_coul;
  double grij,expm2,prefactor,t,erfc;
  double denc;
  int *ilist,*jlist,*numneigh,**firstneigh;

  ecoul = 0.0;
  ev_init(eflag,vflag);

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
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cut_coulsq) {

        r = sqrt(rsq);
        grij = g_ewald * r;
        expm2 = exp(-grij*grij);
        t = 1.0 / (1.0 + EWALD_P*grij);
        erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

        denc = sqrt(lam2[itype][jtype] + rsq);
        prefactor = qqrd2e * lam1[itype][jtype] * qtmp*q[j] / (denc*denc*denc);

        forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;

        fpair = forcecoul;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          prefactor = qqrd2e * lam1[itype][jtype] * qtmp*q[j] / denc;
          ecoul = prefactor*erfc;
          if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             0.0,ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulLongSoft::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(scale,n+1,n+1,"pair:scale");

  memory->create(lambda,n+1,n+1,"pair:lambda");
  memory->create(lam1,n+1,n+1,"pair:lam1");
  memory->create(lam2,n+1,n+1,"pair:lam2");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulLongSoft::settings(int narg, char **arg)
{
  if (narg != 3) error->all(FLERR,"Illegal pair_style command");

  nlambda = utils::numeric(FLERR,arg[0],false,lmp);
  alphac  = utils::numeric(FLERR,arg[1],false,lmp);

  cut_coul = utils::numeric(FLERR,arg[2],false,lmp);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulLongSoft::coeff(int narg, char **arg)
{
  if (narg != 3)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  double lambda_one = utils::numeric(FLERR,arg[2],false,lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      lambda[i][j] = lambda_one;
      scale[i][j] = 1.0;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulLongSoft::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style lj/cut/coul/long requires atom attribute q");

  neighbor->request(this,instance_me);

  cut_coulsq = cut_coul * cut_coul;

  // insure use of KSpace long-range solver, set g_ewald

 if (force->kspace == nullptr)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulLongSoft::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    if (lambda[i][i] != lambda[j][j])
      error->all(FLERR,"Pair coul/cut/soft different lambda values in mix");
    lambda[i][j] = lambda[i][i];
  }

  lam1[i][j] = pow(lambda[i][j], nlambda);
  lam2[i][j] = alphac * (1.0 - lambda[i][j])*(1.0 - lambda[i][j]);

  scale[j][i] = scale[i][j];
  lambda[j][i] = lambda[i][j];
  lam1[j][i] = lam1[i][j];
  lam2[j][i] = lam2[i][j];

  return cut_coul+2.0*qdist;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulLongSoft::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j])
        fwrite(&lambda[i][j],sizeof(double),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulLongSoft::read_restart(FILE *fp)
{
  read_restart_settings(fp);

  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0)
          utils::sfread(FLERR,&lambda[i][j],sizeof(double),1,fp,nullptr,error);
        MPI_Bcast(&lambda[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulLongSoft::write_restart_settings(FILE *fp)
{
  fwrite(&nlambda,sizeof(double),1,fp);
  fwrite(&alphac,sizeof(double),1,fp);

  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulLongSoft::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    utils::sfread(FLERR,&nlambda,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&alphac,sizeof(double),1,fp,nullptr,error);

    utils::sfread(FLERR,&cut_coul,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&offset_flag,sizeof(int),1,fp,nullptr,error);
    utils::sfread(FLERR,&mix_flag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&nlambda,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alphac,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulLongSoft::single(int i, int j, int itype, int jtype,
                            double rsq,
                            double factor_coul, double /*factor_lj*/,
                            double &fforce)
{
  double r,grij,expm2,t,erfc,prefactor;
  double forcecoul,phicoul;
  double denc;

  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    t = 1.0 / (1.0 + EWALD_P*grij);
    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

    denc = sqrt(lam2[itype][jtype] + rsq);
    prefactor = force->qqrd2e * lam1[itype][jtype] * atom->q[i]*atom->q[j] /
      (denc*denc*denc);

    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  } else forcecoul = 0.0;

  fforce = forcecoul;

  if (rsq < cut_coulsq) {
    prefactor = force->qqrd2e * lam1[itype][jtype] * atom->q[i]*atom->q[j] / denc;
    phicoul = prefactor*erfc;
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
  } else phicoul = 0.0;

  return phicoul;
}

/* ---------------------------------------------------------------------- */

void *PairCoulLongSoft::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  dim = 2;
  if (strcmp(str,"scale") == 0) return (void *) scale;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;

  return nullptr;
}
