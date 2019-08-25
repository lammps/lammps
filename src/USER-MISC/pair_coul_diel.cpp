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
   Contributiong authors: Arben Jusufi, Axel Kohlmeyer (Temple U.)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_coul_diel.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulDiel::PairCoulDiel(LAMMPS *lmp) : Pair(lmp) {}

/* ---------------------------------------------------------------------- */

PairCoulDiel::~PairCoulDiel()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(sigmae);
    memory->destroy(rme);
    memory->destroy(offset);
    memory->destroy(cutsq);
    memory->destroy(cut);
    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulDiel::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair;
  double rsq,r,rarg,th,depsdr,epsr,forcecoul,factor_coul;
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

      if (rsq < cutsq[itype][jtype]) {
        r = sqrt(rsq);
        rarg = (r-rme[itype][jtype])/sigmae[itype][jtype];
        th=tanh(rarg);
        epsr=a_eps+b_eps*th;
        depsdr=b_eps * (1.0 - th*th) / sigmae[itype][jtype];

        forcecoul = qqrd2e*qtmp*q[j]*((eps_s*(epsr+r*depsdr)/epsr/epsr) -1.)/rsq;
        fpair = factor_coul*forcecoul/r;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          ecoul = (qqrd2e*qtmp*q[j]*((eps_s/epsr) -1.)/r) - offset[itype][jtype];
          ecoul *= factor_coul;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,0.0,
                             ecoul,fpair,delx,dely,delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairCoulDiel::allocate()
{
  allocated = 1;
  int n = atom->ntypes;

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(sigmae,n+1,n+1,"pair:sigmae");
  memory->create(rme,n+1,n+1,"pair:rme");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulDiel::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

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

void PairCoulDiel::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  eps_s = force->numeric(FLERR,arg[2]);
  double rme_one =force->numeric(FLERR,arg[3]);
  double sigmae_one = force->numeric(FLERR,arg[4]);

  double cut_one = cut_global;
  if (narg == 6) cut_one = force->numeric(FLERR,arg[5]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      sigmae[i][j] = sigmae_one;
      rme[i][j] = rme_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }
  a_eps = 0.5*(5.2+eps_s);
  b_eps = 0.5*(eps_s-5.2);

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulDiel::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/diel requires atom attribute q");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulDiel::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    error->all(FLERR,"for pair style coul/diel, parameters need to be set explicitly for all pairs.");
  }

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  if (offset_flag && (cut[i][j] > 0.0)) {
    double rarg = (cut[i][j]-rme[i][j])/sigmae[i][j];
    double epsr=a_eps+b_eps*tanh(rarg);
    offset[i][j] = qqrd2e*q[i]*q[j]*((eps_s/epsr) -1.)/cut[i][j];
  } else offset[i][j] = 0.0;


  sigmae[j][i] = sigmae[i][j];
  rme[j][i] = rme[i][j];
  offset[j][i] = offset[i][j];
  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDiel::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&rme[i][j],sizeof(double),1,fp);
        fwrite(&sigmae[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDiel::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&rme[i][j],sizeof(double),1,fp);
          fread(&sigmae[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&rme[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigmae[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulDiel::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulDiel::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}

/* ---------------------------------------------------------------------- */

double PairCoulDiel::single(int i, int j, int itype, int jtype,
                           double rsq, double factor_coul, double /*factor_lj*/,
                           double &fforce)
{
  double r, rarg,forcedielec,phidielec;
  double th,epsr,depsdr;
  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

  r=sqrt(rsq);
  rarg = (r-rme[itype][jtype])/sigmae[itype][jtype];
  th = tanh(rarg);
  epsr=a_eps+b_eps*th;
  depsdr=b_eps*(1.-th*th)/sigmae[itype][jtype];

  forcedielec = qqrd2e*q[i]*q[j]*((eps_s*(epsr+r*depsdr)/epsr/epsr) -1.)/rsq;
  fforce = factor_coul*forcedielec/r;

  phidielec = (qqrd2e*q[i]*q[j]*((eps_s/epsr) -1.)/r)- offset[itype][jtype];
  return factor_coul*phidielec;
}
