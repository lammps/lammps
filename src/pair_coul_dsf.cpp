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
   Contributing author: Trung Dac Nguyen (ORNL)
   References: Fennell and Gezelter, JCP 124, 234104 (2006)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_coul_dsf.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_const.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairCoulDSF::PairCoulDSF(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairCoulDSF::~PairCoulDSF()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulDSF::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double r,rsq,r2inv,forcecoul,factor_coul;
  double prefactor,erfcc,erfcd,e_self,t;
  int *ilist,*jlist,*numneigh,**firstneigh;
  
  ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;
  
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
    
    if (evflag) {
      e_self = -(e_shift/2.0 + alpha/MY_PIS) * qtmp*qtmp*qqrd2e;
      ev_tally(i,i,nlocal,0,0.0,e_self,0.0,0.0,0.0,0.0);
    }

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
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
          prefactor = factor_coul * qqrd2e*qtmp*q[j]/r;
          erfcd = exp(-alpha*alpha*rsq);
          t = 1.0 / (1.0 + EWALD_P*alpha*r);
          erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
          forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS * erfcd + 
            r*f_shift) * r;
        }

        fpair = forcecoul * r2inv;
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
          } else ecoul = 0.0;
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

void PairCoulDSF::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulDSF::settings(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Illegal pair_style command");

  alpha = force->numeric(FLERR,arg[0]);
  cut_coul = force->numeric(FLERR,arg[1]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulDSF::coeff(int narg, char **arg)
{
  if (narg != 2) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);
  
  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulDSF::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/dsf requires atom attribute q");

  neighbor->request(this);

  cut_coulsq = cut_coul * cut_coul;
  double erfcc = erfc(alpha*cut_coul); 
  double erfcd = exp(-alpha*alpha*cut_coul*cut_coul);
  f_shift = -(erfcc/cut_coulsq + 2.0/MY_PIS*alpha*erfcd/cut_coul); 
  e_shift = erfcc/cut_coul - f_shift*cut_coul; 
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulDSF::init_one(int i, int j)
{
  return cut_coul;
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDSF::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDSF::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDSF::write_restart_settings(FILE *fp)
{
  fwrite(&alpha,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDSF::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&alpha,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&alpha,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulDSF::single(int i, int j, int itype, int jtype, double rsq,
                           double factor_coul, double factor_lj,
                           double &fforce)
{
  double r2inv,r,erfcc,erfcd,prefactor,t;
  double forcecoul,phicoul;
  
  r2inv = 1.0/rsq;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = factor_coul * force->qqrd2e * atom->q[i]*atom->q[j]/r;
    erfcd = exp(-alpha*alpha*rsq);
    t = 1.0 / (1.0 + EWALD_P*alpha*r);
    erfcc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * erfcd;
	
    forcecoul = prefactor * (erfcc/r + 2.0*alpha/MY_PIS*erfcd + 
      r*f_shift) * r;

    phicoul = prefactor * (erfcc - r*e_shift - rsq*f_shift);
    eng += phicoul;
  } else forcecoul = 0.0;
  
  fforce = forcecoul * r2inv;
      
  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairCoulDSF::extract(const char *str, int &dim)
{
  if (strcmp(str,"cut_coul") == 0) {
    dim = 0;
    return (void *) &cut_coul;
  }
  return NULL;
}
