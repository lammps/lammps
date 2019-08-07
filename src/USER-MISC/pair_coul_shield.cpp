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
   Contributing author: Wengen Ouyang (Tel Aviv University)
   e-mail: w.g.ouyang at gmail dot com

   This is a Coulomb potential described in
   [Maaravi et al, J. Phys. Chem. C 121, 22826-22835 (2017)]
------------------------------------------------------------------------- */

#include "pair_coul_shield.h"
#include <cmath>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "math_special.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairCoulShield::PairCoulShield(LAMMPS *lmp) : Pair(lmp) {
  tap_flag = 1;
}

/* ---------------------------------------------------------------------- */

PairCoulShield::~PairCoulShield()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(sigmae);
    memory->destroy(offset);
    memory->destroy(cutsq);
    memory->destroy(cut);
    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- */

void PairCoulShield::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,ecoul,fpair,Tap,dTap;
  double rsq,r,r3,rarg,th,depsdr,epsr,forcecoul,factor_coul,Vc,fvc;
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

      // only include the interation between different layers
      if (rsq < cutsq[itype][jtype] && atom->molecule[i] != atom->molecule[j]) {
        r = sqrt(rsq);
        r3 = rsq*r;
        rarg = 1.0/sigmae[itype][jtype];
        th = r3 + MathSpecial::cube(rarg);
        epsr = 1.0/pow(th,0.333333333333333333333333);
        depsdr = MathSpecial::square(epsr);
        depsdr *= depsdr;
        Vc = qqrd2e*qtmp*q[j]*epsr;

        // turn on/off taper function
        if (tap_flag) {
          Tap = calc_Tap(r,sqrt(cutsq[itype][jtype]));
          dTap = calc_dTap(r,sqrt(cutsq[itype][jtype]));
        } else {Tap = 1.0; dTap = 0.0;}

        forcecoul = qqrd2e*qtmp*q[j]*r*depsdr;
        fvc = forcecoul*Tap - Vc*dTap/r;
        fpair = factor_coul*fvc;

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }

        if (eflag) {
          if (tap_flag) ecoul = Vc*Tap;
          else ecoul = Vc - offset[itype][jtype];
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

void PairCoulShield::allocate()
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
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairCoulShield::settings(int narg, char **arg)
{
  if (narg < 1 || narg > 2) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);
  if (narg == 2) tap_flag = force->numeric(FLERR,arg[1]);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairCoulShield::coeff(int narg, char **arg)
{
  if (narg < 3 || narg > 4) error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double sigmae_one = force->numeric(FLERR,arg[2]);

  double cut_one = cut_global;
  if (narg == 4) cut_one = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      sigmae[i][j] = sigmae_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients");
}


/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairCoulShield::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,"Pair style coul/shield requires atom attribute q");
  if (!atom->molecule_flag)
    error->all(FLERR,"Pair style coul/shield requires atom attribute molecule");

  neighbor->request(this,instance_me);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairCoulShield::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    error->all(FLERR,"for pair style coul/shield, parameters need to be set explicitly for all pairs.");
  }

  double *q = atom->q;
  double qqrd2e = force->qqrd2e;
  double r,r3,rarg,th,epsr;

  if (offset_flag) {
     r = cut[i][j];
     r3 = r*r*r;
     rarg = 1.0/sigmae[i][j];
     th = r3 + MathSpecial::cube(rarg);
     epsr = 1.0/pow(th,0.333333333333333333);
     offset[i][j] = qqrd2e*q[i]*q[j]*epsr;
  } else offset[i][j] = 0.0;


  sigmae[j][i] = sigmae[i][j];
  offset[j][i] = offset[i][j];
  cut[j][i] = cut[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulShield::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&sigmae[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulShield::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j]) {
        if (me == 0) {
          fread(&sigmae[i][j],sizeof(double),1,fp);
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&sigmae[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairCoulShield::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairCoulShield::read_restart_settings(FILE *fp)
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

double PairCoulShield::single(int i, int j, int itype, int jtype,
                           double rsq, double factor_coul, double /*factor_lj*/,
                           double &fforce)
{
  double r, rarg,Vc,fvc,forcecoul,phishieldec;
  double r3,th,epsr,depsdr,Tap,dTap;
  double *q = atom->q;
  double qqrd2e = force->qqrd2e;

   r = sqrt(rsq);
   r3 = rsq*r;
   rarg = 1.0/sigmae[itype][jtype];
   th = r3 + MathSpecial::cube(rarg);
   epsr = 1.0/pow(th,0.333333333333333333);
   depsdr = epsr*epsr;
   depsdr *= depsdr;
   Vc = qqrd2e*q[i]*q[j]*epsr;

   // turn on/off taper function
   if (tap_flag) {
     Tap = calc_Tap(r,sqrt(cutsq[itype][jtype]));
     dTap = calc_dTap(r,sqrt(cutsq[itype][jtype]));
   } else {Tap = 1.0; dTap = 0.0;}

   forcecoul = qqrd2e*q[i]*q[j]*r*depsdr;
   fvc = forcecoul*Tap - Vc*dTap/r;
   fforce = factor_coul*fvc;

  if (tap_flag) phishieldec = factor_coul*Vc*Tap;
  else phishieldec = Vc - offset[itype][jtype];
  return factor_coul*phishieldec;
}
