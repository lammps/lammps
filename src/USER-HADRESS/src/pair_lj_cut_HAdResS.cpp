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
   Contributing author: Maziar Heidari,
                        Robinson Cortes-Huerto,
                        Davide Donadio and
                        Raffaello Potestio
                        (Max Planck Institute for Polymer Research, 2015-2016)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut_HAdResS.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

PairLJCutHADRESS::PairLJCutHADRESS(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJCutHADRESS::~PairLJCutHADRESS()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lj1);
    memory->destroy(lj2);
    memory->destroy(lj3);
    memory->destroy(lj4);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutHADRESS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj, forcelj11, forcelj22,factor_lj,gradLambdai,gradLambdaj,V11,V22;
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

  /* Start of MaHe Contribution */
  double x0lo = domain->boxlo[0];
  double x0hi = domain->boxhi[0];
  double x0BoxSize = x0hi - x0lo;

  double d1 = x0lo + SizeRegionA/2.0;
  double d2 = x0lo + x0BoxSize - (SizeRegionA/2.0 + dHyb);
  double xtmpj, iLambda, jLambda, ijLambda;

  /* End of MaHe Contribution */


  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    if (!THIN){
/* Start of MaHe Contribution */

    if(xtmp <= d1)iLambda = 1.0;
    else if(xtmp >= d1 && xtmp < d1+dHyb) iLambda = cos(MY_PI*(xtmp-d1)/(2.0*dHyb)) * cos(MY_PI*(xtmp-d1)/(2.0*dHyb));
    else if(xtmp >=d1+dHyb && xtmp < d2)iLambda = 0.0;
    else if(xtmp >= d2 && xtmp < d2+dHyb) iLambda = cos(MY_PI*(xtmp-d2+dHyb)/(2.0*dHyb)) * cos(MY_PI*(xtmp-d2+dHyb)/(2.0*dHyb));
    else if(xtmp >=d2+dHyb)iLambda = 1.0;

    if(xtmp <= d1)gradLambdai = 0.0;
    else if(xtmp >= d1 && xtmp < d1+dHyb) gradLambdai = -(MY_PI/dHyb)*sin(MY_PI*(xtmp-d1)/(2.0*dHyb)) * cos(MY_PI*(xtmp-d1)/(2.0*dHyb));
    else if(xtmp >=d1+dHyb && xtmp < d2)gradLambdai = 0.0;
    else if(xtmp >= d2 && xtmp < d2+dHyb) gradLambdai = -(MY_PI/dHyb)*sin(MY_PI*(xtmp-d2+dHyb)/(2.0*dHyb)) * cos(MY_PI*(xtmp-d2+dHyb)/(2.0*dHyb));
    else if(xtmp >=d2+dHyb)gradLambdai = 0.0;

/* End of MaHe Contribution */

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];

      xtmpj = x[j][0];

      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        /* Start of MaHe Contribution */
        if(xtmpj <= d1)jLambda = 1.0;
        else if(xtmpj >= d1 && xtmpj < d1+dHyb) jLambda = cos(MY_PI*(xtmpj-d1)/(2.0*dHyb)) * cos(MY_PI*(xtmpj-d1)/(2.0*dHyb));
        else if(xtmpj >=d1+dHyb && xtmpj < d2)jLambda = 0.0;
        else if(xtmpj >= d2 && xtmpj < d2+dHyb) jLambda = cos(MY_PI*(xtmpj-d2+dHyb)/(2.0*dHyb)) * cos(MY_PI*(xtmpj-d2+dHyb)/(2.0*dHyb));
        else if(xtmpj >=d2+dHyb)jLambda = 1.0;
        /* End of MaHe Contribution */

        if((iLambda==0 && jLambda==0 && lj1[2][2]!=0) || (iLambda==1 && jLambda==1 && lj1[1][1]!=0)){

            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            if(jLambda==1)forcelj = r6inv * (lj1[1][1]*r6inv - lj2[1][1]);
            if(jLambda==0)forcelj = r6inv * (lj1[2][2]*r6inv - lj2[2][2]);

            fpair = factor_lj*forcelj*r2inv;


            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;
            /* End of MaHe Contribution */


            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*fpair;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;
            }

            if (eflag) {
                if(jLambda==1)evdwl = r6inv*(lj3[1][1]*r6inv-lj4[1][1]) -
                        offset[1][1];
                if(jLambda==0)evdwl = r6inv*(lj3[2][2]*r6inv-lj4[2][2]) -
                        offset[2][2];
              evdwl *= factor_lj;
            }

            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fpair,delx,dely,delz);

        }
        else{


            if(xtmpj <= d1)gradLambdaj = 0.0;
            else if(xtmpj >= d1 && xtmpj < d1+dHyb) gradLambdaj = -(MY_PI/dHyb)*sin(MY_PI*(xtmpj-d1)/(2.0*dHyb)) * cos(MY_PI*(xtmpj-d1)/(2.0*dHyb));
            else if(xtmpj >=d1+dHyb && xtmpj < d2)gradLambdaj = 0.0;
            else if(xtmpj >= d2 && xtmpj < d2+dHyb) gradLambdaj = -(MY_PI/dHyb)*sin(MY_PI*(xtmpj-d2+dHyb)/(2.0*dHyb)) * cos(MY_PI*(xtmpj-d2+dHyb)/(2.0*dHyb));
            else if(xtmpj >=d2+dHyb)gradLambdaj = 0.0;

            ijLambda = 0.5 * (iLambda + jLambda);
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj11 = r6inv * (lj1[1][1]*r6inv - lj2[1][1]);
            forcelj22 = r6inv * (lj1[2][2]*r6inv - lj2[2][2]);

            fpair = factor_lj*(forcelj11*ijLambda+forcelj22*(1-ijLambda))*r2inv;

            V11 = r6inv*(lj3[1][1]*r6inv-lj4[1][1]);
            V22 = r6inv*(lj3[2][2]*r6inv-lj4[2][2]);

            f[i][0] += delx*fpair - 0.5*(V11-V22)*gradLambdai;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;
            /* End of MaHe Contribution */


            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*fpair + 0.5*(V11-V22)*gradLambdaj;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;
            }

            if (eflag) {
                evdwl = (ijLambda*V11 + (1-ijLambda)*V22);
              evdwl *= factor_lj;
            }

            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fpair,delx,dely,delz);

        }

        }
    }
    }
    else{
        for (jj = 0; jj < jnum; jj++) {
          j = jlist[jj];
          factor_lj = special_lj[sbmask(j)];
          j &= NEIGHMASK;

          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];

          xtmpj = x[j][0];

          rsq = delx*delx + dely*dely + delz*delz;
          jtype = type[j];

          if (rsq < cutsq[itype][jtype]) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj11 = r6inv * (lj1[1][1]*r6inv - lj2[1][1]);
            forcelj22 = r6inv * (lj1[2][2]*r6inv - lj2[2][2]);

            lambda_THIN += lmp->update->atimestep * lambda_Increment;
            if(lambda_THIN > 1.0) lambda_THIN = 1.0;

            fpair = factor_lj*(forcelj11*lambda_THIN+forcelj22*(1-lambda_THIN))*r2inv;

            /* Start of MaHe Contribution */
            f[i][0] += delx*fpair;
            f[i][1] += dely*fpair;
            f[i][2] += delz*fpair;
            /* End of MaHe Contribution */


            if (newton_pair || j < nlocal) {
              f[j][0] -= delx*fpair;
              f[j][1] -= dely*fpair;
              f[j][2] -= delz*fpair;
            }

            if (eflag) {
                evdwl =  r6inv*(lambda_THIN*(lj3[1][1]*r6inv-lj4[2][2]) + (1-lambda_THIN)*(lj3[1][1]*r6inv-lj4[2][2])) -
              offset[itype][jtype];
              evdwl *= factor_lj;
            }

//            if (eflag) {
              //Free_Energy_Diff =  r6inv*(lj3AA*r6inv-lj4AA-lj3BB*r6inv-lj4BB);
//            }


            if (evflag) ev_tally(i,j,nlocal,newton_pair,
                                 evdwl,0.0,fpair,delx,dely,delz);
          }
        }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLJCutHADRESS::compute_inner()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj,rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = listinner->inum;
  ilist = listinner->ilist;
  numneigh = listinner->numneigh;
  firstneigh = listinner->firstneigh;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];

  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

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

      if (rsq < cut_out_off_sq) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        jtype = type[j];
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair *= 1.0 - rsw*rsw*(3.0 - 2.0*rsw);
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutHADRESS::compute_middle()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj,rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = listmiddle->inum;
  ilist = listmiddle->ilist;
  numneigh = listmiddle->numneigh;
  firstneigh = listmiddle->firstneigh;

  double cut_in_off = cut_respa[0];
  double cut_in_on = cut_respa[1];
  double cut_out_on = cut_respa[2];
  double cut_out_off = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

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

      if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
        r2inv = 1.0/rsq;
        r6inv = r2inv*r2inv*r2inv;
        jtype = type[j];
        forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
        fpair = factor_lj*forcelj*r2inv;
        if (rsq < cut_in_on_sq) {
          rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
          fpair *= rsw*rsw*(3.0 - 2.0*rsw);
        }
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair *= 1.0 + rsw*rsw*(2.0*rsw - 3.0);
        }

        f[i][0] += delx*fpair;
        f[i][1] += dely*fpair;
        f[i][2] += delz*fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx*fpair;
          f[j][1] -= dely*fpair;
          f[j][2] -= delz*fpair;
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutHADRESS::compute_outer(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,r2inv,r6inv,forcelj,factor_lj,rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = listouter->inum;
  ilist = listouter->ilist;
  numneigh = listouter->numneigh;
  firstneigh = listouter->firstneigh;

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

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
        if (rsq > cut_in_off_sq) {
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
          fpair = factor_lj*forcelj*r2inv;
          if (rsq < cut_in_on_sq) {
            rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
            fpair *= rsw*rsw*(3.0 - 2.0*rsw);
          }

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
        }

        if (eflag) {
          r2inv = 1.0/rsq;
          r6inv = r2inv*r2inv*r2inv;
          evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
            offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (vflag) {
          if (rsq <= cut_in_off_sq) {
            r2inv = 1.0/rsq;
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
            fpair = factor_lj*forcelj*r2inv;
          } else if (rsq < cut_in_on_sq)
            fpair = factor_lj*forcelj*r2inv;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,0.0,fpair,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutHADRESS::allocate()
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
  memory->create(lj1,n+1,n+1,"pair:lj1");
  memory->create(lj2,n+1,n+1,"pair:lj2");
  memory->create(lj3,n+1,n+1,"pair:lj3");
  memory->create(lj4,n+1,n+1,"pair:lj4");
  memory->create(offset,n+1,n+1,"pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutHADRESS::settings(int narg, char **arg)
{
  //if (narg != 5) error->all(FLERR,"Illegal pair_style command");

  cut_global = force->numeric(FLERR,arg[0]);

  THIN  = force->numeric(FLERR,arg[1]);

  dHyb  = force->numeric(FLERR,arg[2]);
  SizeRegionA  = force->numeric(FLERR,arg[3]);

  lambda_Increment  = force->numeric(FLERR,arg[4]);
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

void PairLJCutHADRESS::coeff(int narg, char **arg)
{
    //for(int i=0;i<narg;i++) printf("arg %s\n",arg[i]);
  if (narg < 4 || narg > 5)
    error->all(FLERR,"Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);


  double cut_one = cut_global;
  if (narg == 5) cut_one = force->numeric(FLERR,arg[4]);


  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
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

void PairLJCutHADRESS::init_style()
{
  // request regular or rRESPA neighbor lists

  int irequest;

  if (update->whichflag == 1 && strstr(update->integrate_style,"respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;

    if (respa == 0) irequest = neighbor->request(this,instance_me);
    else if (respa == 1) {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    } else {
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 1;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respainner = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 2;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respamiddle = 1;
      irequest = neighbor->request(this,instance_me);
      neighbor->requests[irequest]->id = 3;
      neighbor->requests[irequest]->half = 0;
      neighbor->requests[irequest]->respaouter = 1;
    }

  } else irequest = neighbor->request(this,instance_me);

  // set rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0)
    cut_respa = ((Respa *) update->integrate)->cutoff;
  else cut_respa = NULL;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   regular or rRESPA
------------------------------------------------------------------------- */

void PairLJCutHADRESS::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listinner = ptr;
  else if (id == 2) listmiddle = ptr;
  else if (id == 3) listouter = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCutHADRESS::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  //printf("i=%d, j=%d, sigma= %f\n",i,j,sigma[i][j]);

  lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j],6.0);
  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],12.0);
  lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j],6.0);


  if (offset_flag) {
    double ratio = sigma[i][j] / cut[i][j];
    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio,12.0) - pow(ratio,6.0));
  } else offset[i][j] = 0.0;

  lj1[j][i] = lj1[i][j];
  lj2[j][i] = lj2[i][j];
  lj3[j][i] = lj3[i][j];
  lj4[j][i] = lj4[i][j];
  offset[j][i] = offset[i][j];

  // check interior rRESPA cutoff

  if (cut_respa && cut[i][j] < cut_respa[3])
    error->all(FLERR,"Pair cutoff < Respa interior cutoff");

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

    double sig2 = sigma[i][j]*sigma[i][j];
    double sig6 = sig2*sig2*sig2;
    double rc3 = cut[i][j]*cut[i][j]*cut[i][j];
    double rc6 = rc3*rc3;
    double rc9 = rc3*rc6;
    etail_ij = 8.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (sig6 - 3.0*rc6) / (9.0*rc9);
    ptail_ij = 16.0*MY_PI*all[0]*all[1]*epsilon[i][j] *
      sig6 * (2.0*sig6 - 3.0*rc6) / (9.0*rc9);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutHADRESS::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&epsilon[i][j],sizeof(double),1,fp);
        fwrite(&sigma[i][j],sizeof(double),1,fp);
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutHADRESS::read_restart(FILE *fp)
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
          fread(&cut[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutHADRESS::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
  fwrite(&tail_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutHADRESS::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    fread(&cut_global,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
    fread(&tail_flag,sizeof(int),1,fp);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
  MPI_Bcast(&tail_flag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCutHADRESS::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g\n",i,epsilon[i][i],sigma[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCutHADRESS::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJCutHADRESS::single(int i, int j, int itype, int jtype, double rsq,
                         double factor_coul, double factor_lj,
                         double &fforce)
{
  double r2inv,r6inv,forcelj,philj;

  r2inv = 1.0/rsq;
  r6inv = r2inv*r2inv*r2inv;
  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
  fforce = factor_lj*forcelj*r2inv;

  philj = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
    offset[itype][jtype];
  return factor_lj*philj;
}

/* ---------------------------------------------------------------------- */

void *PairLJCutHADRESS::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  return NULL;
}
