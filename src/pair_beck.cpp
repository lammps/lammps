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
   Contributing author: Jonathan Zimmerman (Sandia)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "pair_beck.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "math_special.h"

using namespace LAMMPS_NS;
using namespace MathSpecial;

/* ---------------------------------------------------------------------- */

PairBeck::PairBeck(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairBeck::~PairBeck()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(AA);
    memory->destroy(BB);
    memory->destroy(aa);
    memory->destroy(alpha);
    memory->destroy(beta);
  }
}

/* ---------------------------------------------------------------------- */

void PairBeck::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r5,force_beck,factor_lj;
  double r,rinv;
  double aaij,alphaij,betaij;
  double term1,term1inv,term2,term3,term4,term5,term6;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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
        r = sqrt(rsq);
        r5 = rsq*rsq*r;
        aaij = aa[itype][jtype];
        alphaij = alpha[itype][jtype];
        betaij = beta[itype][jtype];
        term1 = aaij*aaij + rsq;
        term2 = powint(term1,-5);
        term3 = 21.672 + 30.0*aaij*aaij + 6.0*rsq;
        term4 = alphaij + r5*betaij;
        term5 = alphaij + 6.0*r5*betaij;
        rinv  = 1.0/r;
        force_beck = AA[itype][jtype]*exp(-1.0*r*term4)*term5;
        force_beck -= BB[itype][jtype]*r*term2*term3;

        fpair = factor_lj*force_beck*rinv;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          term6 = powint(term1,-3);
          term1inv = 1.0/term1;
          evdwl = AA[itype][jtype]*exp(-1.0*r*term4);
          evdwl -= BB[itype][jtype]*term6*(1.0+(2.709+3.0*aaij*aaij)*term1inv);
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

void PairBeck::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");

  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(AA,n+1,n+1,"pair:AA");
  memory->create(BB,n+1,n+1,"pair:BB");
  memory->create(aa,n+1,n+1,"pair:aa");
  memory->create(alpha,n+1,n+1,"pair:alpha");
  memory->create(beta,n+1,n+1,"pair:beta");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairBeck::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut[i][j] = cut_global;
        }
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairBeck::coeff(int narg, char **arg)
{
  if (narg != 7 && narg != 8)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double AA_one = force->numeric(FLERR,arg[2]);
  double BB_one = force->numeric(FLERR,arg[3]);
  double aa_one = force->numeric(FLERR,arg[4]);
  double alpha_one = force->numeric(FLERR,arg[5]);
  double beta_one = force->numeric(FLERR,arg[6]);

  double cut_one = cut_global;
  if (narg == 8) cut_one = force->numeric(FLERR,arg[7]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      AA[i][j] = AA_one;
      BB[i][j] = BB_one;
      aa[i][j] = aa_one;
      alpha[i][j] = alpha_one;
      beta[i][j] = beta_one;
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

double PairBeck::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  AA[j][i] = AA[i][j];
  BB[j][i] = BB[i][j];
  aa[j][i] = aa[i][j];
  alpha[j][i] = alpha[i][j];
  beta[j][i] = beta[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBeck::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&AA[i][j],sizeof(double),1,fp);
        fwrite(&BB[i][j],sizeof(double),1,fp);
        fwrite(&aa[i][j],sizeof(double),1,fp);
        fwrite(&alpha[i][j],sizeof(double),1,fp);
        fwrite(&beta[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBeck::read_restart(FILE *fp)
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
          fread(&AA[i][j],sizeof(double),1,fp);
          fread(&BB[i][j],sizeof(double),1,fp);
          fread(&aa[i][j],sizeof(double),1,fp);
          fread(&alpha[i][j],sizeof(double),1,fp);
          fread(&beta[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&AA[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&BB[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&aa[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&alpha[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&beta[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairBeck::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairBeck::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairBeck::single(int i, int j, int itype, int jtype,
                                  double rsq,
                                  double factor_coul, double factor_lj,
                                  double &fforce)
{
  double phi_beck,r,rinv;
  double r5,force_beck;
  double aaij,alphaij,betaij;
  double term1,term1inv,term2,term3,term4,term5,term6;

  r = sqrt(rsq);
  r5 = rsq*rsq*r;
  aaij = aa[itype][jtype];
  alphaij = alpha[itype][jtype];
  betaij = beta[itype][jtype];
  term1 = aaij*aaij + rsq;
  term2 = powint(term1,-5);
  term3 = 21.672 + 30.0*aaij*aaij + 6.0*rsq;
  term4 = alphaij + r5*betaij;
  term5 = alphaij + 6.0*r5*betaij;
  rinv  = 1.0/r;
  force_beck = AA[itype][jtype]*exp(-1.0*r*term4)*term5;
  force_beck -= BB[itype][jtype]*r*term2*term3;
  fforce = factor_lj*force_beck*rinv;

  term6 = powint(term1,-3);
  term1inv = 1.0/term1;
  phi_beck = AA[itype][jtype]*exp(-1.0*r*term4);
  phi_beck -= BB[itype][jtype]*term6*(1.0+(2.709+3.0*aaij*aaij)*term1inv);

  return factor_lj*phi_beck;
}
