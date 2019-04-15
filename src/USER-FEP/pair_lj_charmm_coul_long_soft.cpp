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
   Soft-core version: Agilio Padua (Univ Blaise Pascal & CNRS)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_lj_charmm_coul_long_soft.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
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

PairLJCharmmCoulLongSoft::PairLJCharmmCoulLongSoft(LAMMPS *lmp) : Pair(lmp)
{
  respa_enable = 1;
  ewaldflag = pppmflag = 1;
  implicit = 0;
  mix_flag = ARITHMETIC;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLongSoft::~PairLJCharmmCoulLongSoft()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(epsilon);
    memory->destroy(sigma);
    memory->destroy(lambda);
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

void PairLJCharmmCoulLongSoft::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,forcecoul,forcelj,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double philj,switch1,switch2;
  double denc, denlj, r4sig6;
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
      jtype = type[j];

      if (rsq < cut_bothsq) {

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          denc = sqrt(lj4[itype][jtype] + rsq);
          prefactor = qqrd2e * lj1[itype][jtype] * qtmp*q[j] / (denc*denc*denc);

          forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq) {
          r4sig6 = rsq*rsq / lj2[itype][jtype];
          denlj = lj3[itype][jtype] + rsq*r4sig6;
          forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
            (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));

          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
              (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
            switch2 = 12.0 * (cut_ljsq-rsq) * (rsq-cut_lj_innersq) / denom_lj;
            philj = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
              (1.0/(denlj*denlj) - 1.0/denlj);
            forcelj = forcelj*switch1 + philj*switch2;
          }
        } else forcelj = 0.0;

        fpair = forcecoul + factor_lj*forcelj;

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
            prefactor = qqrd2e * lj1[itype][jtype] * qtmp*q[j] / denc;
            ecoul = prefactor*erfc;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;

          if (rsq < cut_ljsq) {
            evdwl = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
              (1.0/(denlj*denlj) - 1.0/denlj);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
              evdwl *= switch1;
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

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::compute_inner()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,forcecoul,forcelj,factor_coul,factor_lj;
  double rsw;
  double denc, denlj, r4sig6;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum_inner;
  ilist = list->ilist_inner;
  numneigh = list->numneigh_inner;
  firstneigh = list->firstneigh_inner;

  double cut_out_on = cut_respa[0];
  double cut_out_off = cut_respa[1];

  double cut_out_diff = cut_out_off - cut_out_on;
  double cut_out_on_sq = cut_out_on*cut_out_on;
  double cut_out_off_sq = cut_out_off*cut_out_off;

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

      if (rsq < cut_out_off_sq) {
        jtype = type[j];

        denc = sqrt(lj4[itype][jtype] + rsq);
        forcecoul = qqrd2e * lj1[itype][jtype] * qtmp*q[j] / (denc*denc*denc);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        r4sig6 = rsq*rsq / lj2[itype][jtype];
        denlj = lj3[itype][jtype] + rsq*r4sig6;
        forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
          (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));

        fpair = forcecoul + factor_lj*forcelj;

        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
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

void PairLJCharmmCoulLongSoft::compute_middle()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,forcecoul,forcelj,factor_coul,factor_lj;
  double philj,switch1,switch2;
  double rsw;
  double denc, denlj, r4sig6;
  int *ilist,*jlist,*numneigh,**firstneigh;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum_middle;
  ilist = list->ilist_middle;
  numneigh = list->numneigh_middle;
  firstneigh = list->firstneigh_middle;

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

      if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
        jtype = type[j];

        denc = sqrt(lj4[itype][jtype] + rsq);
        forcecoul = qqrd2e * lj1[itype][jtype] * qtmp*q[j] / (denc*denc*denc);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        r4sig6 = rsq*rsq / lj2[itype][jtype];
        denlj = lj3[itype][jtype] + rsq*r4sig6;
        forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
          (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));

        if (rsq > cut_lj_innersq) {
          switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
            (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
          switch2 = 12.0 * (cut_ljsq-rsq) * (rsq-cut_lj_innersq) / denom_lj;
          philj = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
              (1.0/(denlj*denlj) - 1.0/denlj);
          forcelj = forcelj*switch1 + philj*switch2;
        }

        fpair = forcecoul + factor_lj*forcelj;

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

void PairLJCharmmCoulLongSoft::compute_outer(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double r,rsq,forcecoul,forcelj,factor_coul,factor_lj;
  double grij,expm2,fprefactor,eprefactor,t,erfc;
  double philj,switch1,switch2;
  double rsw;
  double denc, denlj, r4sig6;
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

  double cut_in_off = cut_respa[2];
  double cut_in_on = cut_respa[3];

  double cut_in_diff = cut_in_on - cut_in_off;
  double cut_in_off_sq = cut_in_off*cut_in_off;
  double cut_in_on_sq = cut_in_on*cut_in_on;

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
        jtype = type[j];

        if (rsq < cut_coulsq) {
          r = sqrt(rsq);
          grij = g_ewald * r;
          expm2 = exp(-grij*grij);
          t = 1.0 / (1.0 + EWALD_P*grij);
          erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          denc = sqrt(lj4[itype][jtype] + rsq);
          fprefactor = qqrd2e * lj1[itype][jtype] * qtmp*q[j] /
            (denc*denc*denc);

          forcecoul = fprefactor * (erfc + EWALD_F*grij*expm2 - 1.0);

          if (rsq > cut_in_off_sq) {
            if (rsq < cut_in_on_sq) {
              rsw = (r - cut_in_off)/cut_in_diff;
              forcecoul += fprefactor*rsw*rsw*(3.0 - 2.0*rsw);
              if (factor_coul < 1.0)
                forcecoul -=
                  (1.0-factor_coul)*fprefactor*rsw*rsw*(3.0 - 2.0*rsw);
            } else {
              forcecoul += fprefactor;
              if (factor_coul < 1.0)
                forcecoul -= (1.0-factor_coul)*fprefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq && rsq > cut_in_off_sq) {
          r4sig6 = rsq*rsq / lj2[itype][jtype];
          denlj = lj3[itype][jtype] + rsq*r4sig6;
          forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
            (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));

          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
              (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
            switch2 = 12.0 * (cut_ljsq-rsq) * (rsq-cut_lj_innersq) / denom_lj;
            philj = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
              (1.0/(denlj*denlj) - 1.0/denlj);
            forcelj = forcelj*switch1 + philj*switch2;
          }
          if (rsq < cut_in_on_sq) {
            rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
            forcelj *= rsw*rsw*(3.0 - 2.0*rsw);
          }
        } else forcelj = 0.0;

        fpair = forcecoul + forcelj;

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
            eprefactor = qqrd2e * lj1[itype][jtype] * qtmp*q[j] / denc;
            ecoul = eprefactor*erfc;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*eprefactor;
          } else ecoul = 0.0;

          if (rsq < cut_ljsq) {
            r4sig6 = rsq*rsq / lj2[itype][jtype];
            denlj = lj3[itype][jtype] + rsq*r4sig6;
            evdwl = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
              (1.0/(denlj*denlj) - 1.0/denlj);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
              evdwl *= switch1;
            }
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (vflag) {
          if (rsq < cut_coulsq) {
            forcecoul = fprefactor * (erfc + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*fprefactor;
          } else forcecoul = 0.0;

          if (rsq <= cut_in_off_sq) {
            r4sig6 = rsq*rsq / lj2[itype][jtype];
            denlj = lj3[itype][jtype] + rsq*r4sig6;
            forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
              (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));

            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
              switch2 = 12.0 * (cut_ljsq-rsq) * (rsq-cut_lj_innersq) / denom_lj;
              philj = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
                (1.0/(denlj*denlj) - 1.0/denlj);
              forcelj = forcelj*switch1 + philj*switch2;
            }
          } else if (rsq <= cut_in_on_sq) {
            r4sig6 = rsq*rsq / lj2[itype][jtype];
            denlj = lj3[itype][jtype] + rsq*r4sig6;
            forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
              (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
              switch2 = 12.0 * (cut_ljsq-rsq) * (rsq-cut_lj_innersq) / denom_lj;
              philj = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
                (1.0/(denlj*denlj) - 1.0/denlj);
              forcelj = forcelj*switch1 + philj*switch2;
            }
          }
          fpair = forcecoul + factor_lj*forcelj;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::allocate()
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
  memory->create(lambda,n+1,n+1,"pair:lambda");
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

void PairLJCharmmCoulLongSoft::settings(int narg, char **arg)
{
  if (narg != 5 && narg != 6) error->all(FLERR,"Illegal pair_style command");

  nlambda = force->numeric(FLERR,arg[0]);
  alphalj = force->numeric(FLERR,arg[1]);
  alphac  = force->numeric(FLERR,arg[2]);

  cut_lj_inner = force->numeric(FLERR,arg[3]);
  cut_lj = force->numeric(FLERR,arg[4]);
  if (narg == 5) cut_coul = cut_lj;
  else cut_coul = force->numeric(FLERR,arg[5]);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::coeff(int narg, char **arg)
{
  if (narg != 5 && narg != 7) error->all(FLERR,"Illegal pair_coeff command");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(FLERR,arg[0],atom->ntypes,ilo,ihi);
  force->bounds(FLERR,arg[1],atom->ntypes,jlo,jhi);

  double epsilon_one = force->numeric(FLERR,arg[2]);
  double sigma_one = force->numeric(FLERR,arg[3]);
  double lambda_one = force->numeric(FLERR,arg[4]);

  double eps14_one = epsilon_one;
  double sigma14_one = sigma_one;
  if (narg == 7) {
    eps14_one = force->numeric(FLERR,arg[5]);
    sigma14_one = force->numeric(FLERR,arg[6]);
  }

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      epsilon[i][j] = epsilon_one;
      sigma[i][j] = sigma_one;
      lambda[i][j] = lambda_one;
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

void PairLJCharmmCoulLongSoft::init_style()
{
  if (!atom->q_flag)
    error->all(FLERR,
               "Pair style lj/charmm/coul/long/soft requires atom attribute q");

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

  // require cut_lj_inner < cut_lj

  if (cut_lj_inner >= cut_lj)
    error->all(FLERR,"Pair inner cutoff >= Pair outer cutoff");

  cut_lj_innersq = cut_lj_inner * cut_lj_inner;
  cut_ljsq = cut_lj * cut_lj;
  cut_coulsq = cut_coul * cut_coul;
  cut_bothsq = MAX(cut_ljsq,cut_coulsq);

  denom_lj = (cut_ljsq-cut_lj_innersq) * (cut_ljsq-cut_lj_innersq) *
    (cut_ljsq-cut_lj_innersq);

  // set & error check interior rRESPA cutoffs

  if (strstr(update->integrate_style,"respa") &&
      ((Respa *) update->integrate)->level_inner >= 0) {
    cut_respa = ((Respa *) update->integrate)->cutoff;
    if (MIN(cut_lj,cut_coul) < cut_respa[3])
      error->all(FLERR,"Pair cutoff < Respa interior cutoff");
    if (cut_lj_inner < cut_respa[1])
      error->all(FLERR,"Pair inner cutoff < Respa interior cutoff");
  } else cut_respa = NULL;

  // insure use of KSpace long-range solver, set g_ewald

  if (force->kspace == NULL)
    error->all(FLERR,"Pair style requires a KSpace style");
  g_ewald = force->kspace->g_ewald;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLJCharmmCoulLongSoft::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    epsilon[i][j] = mix_energy(epsilon[i][i],epsilon[j][j],
                               sigma[i][i],sigma[j][j]);
    sigma[i][j] = mix_distance(sigma[i][i],sigma[j][j]);
    if (lambda[i][i] != lambda[j][j])
      error->all(FLERR,"Pair lj/charmm/coul/long/soft different lambda values in mix");
    lambda[i][j] = lambda[i][i];
    eps14[i][j] = mix_energy(eps14[i][i],eps14[j][j],
                               sigma14[i][i],sigma14[j][j]);
    sigma14[i][j] = mix_distance(sigma14[i][i],sigma14[j][j]);
  }

  double cut = MAX(cut_lj,cut_coul);

  lj1[i][j] = pow(lambda[i][j], nlambda);
  lj2[i][j] = pow(sigma[i][j], 6.0);
  lj3[i][j] = alphalj * (1.0 - lambda[i][j])*(1.0 - lambda[i][j]);
  lj4[i][j] = alphac  * (1.0 - lambda[i][j])*(1.0 - lambda[i][j]);

  // 1-4 interactions unaffected (they're part of the dihedral term)
  lj14_1[i][j] = 48.0 * eps14[i][j] * pow(sigma14[i][j],12.0);
  lj14_2[i][j] = 24.0 * eps14[i][j] * pow(sigma14[i][j],6.0);
  lj14_3[i][j] = 4.0 * eps14[i][j] * pow(sigma14[i][j],12.0);
  lj14_4[i][j] = 4.0 * eps14[i][j] * pow(sigma14[i][j],6.0);

  epsilon[j][i] = epsilon[i][j];
  sigma[j][i] = sigma[i][j];
  lambda[j][i] = lambda[i][j];
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
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::write_restart(FILE *fp)
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
        fwrite(&eps14[i][j],sizeof(double),1,fp);
        fwrite(&sigma14[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::read_restart(FILE *fp)
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
          fread(&lambda[i][j],sizeof(double),1,fp);
          fread(&eps14[i][j],sizeof(double),1,fp);
          fread(&sigma14[i][j],sizeof(double),1,fp);
        }
        MPI_Bcast(&epsilon[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&lambda[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&eps14[i][j],1,MPI_DOUBLE,0,world);
        MPI_Bcast(&sigma14[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::write_restart_settings(FILE *fp)
{
  fwrite(&nlambda,sizeof(double),1,fp);
  fwrite(&alphalj,sizeof(double),1,fp);
  fwrite(&alphac,sizeof(double),1,fp);

  fwrite(&cut_lj_inner,sizeof(double),1,fp);
  fwrite(&cut_lj,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&nlambda,sizeof(double),1,fp);
    fread(&alphalj,sizeof(double),1,fp);
    fread(&alphac,sizeof(double),1,fp);

    fread(&cut_lj_inner,sizeof(double),1,fp);
    fread(&cut_lj,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }

  MPI_Bcast(&nlambda,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alphalj,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&alphac,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&cut_lj_inner,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_lj,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);
}


/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d %g %g %g %g %g\n",i,epsilon[i][i],sigma[i][i],
            lambda[i][i],eps14[i][i],sigma14[i][i]);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairLJCharmmCoulLongSoft::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g %g %g %g %g\n",i,j,epsilon[i][j],sigma[i][j],
              lambda[i][j],eps14[i][j],sigma14[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairLJCharmmCoulLongSoft::single(int i, int j, int itype, int jtype,
                                        double rsq,
                                        double factor_coul, double factor_lj,
                                        double &fforce)
{
  double r,grij,expm2,t,erfc,prefactor;
  double switch1,switch2,forcecoul,forcelj,phicoul,philj;
  double denc, denlj, r4sig6;

  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    grij = g_ewald * r;
    expm2 = exp(-grij*grij);
    t = 1.0 / (1.0 + EWALD_P*grij);
    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

    denc = sqrt(lj4[itype][jtype] + rsq);
    prefactor = force->qqrd2e * lj1[itype][jtype] * atom->q[i]*atom->q[j] /
      (denc*denc*denc);

    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  } else forcecoul = 0.0;

  if (rsq < cut_ljsq) {
    r4sig6 = rsq*rsq / lj2[itype][jtype];
    denlj = lj3[itype][jtype] + rsq*r4sig6;
    forcelj = lj1[itype][jtype] * epsilon[itype][jtype] *
      (48.0*r4sig6/(denlj*denlj*denlj) - 24.0*r4sig6/(denlj*denlj));
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
        (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      switch2 = 12.0 * (cut_ljsq-rsq) * (rsq-cut_lj_innersq) / denom_lj;
      philj = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
        (1.0/(denlj*denlj) - 1.0/denlj);
      forcelj = forcelj*switch1 + philj*switch2;
    }
  } else forcelj = 0.0;
  fforce = forcecoul + factor_lj*forcelj;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    prefactor = force->qqrd2e * lj1[itype][jtype] * atom->q[i]*atom->q[j] / denc;
    phicoul = prefactor*erfc;
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    eng += phicoul;
  }

  if (rsq < cut_ljsq) {
    philj = lj1[itype][jtype] * 4.0 * epsilon[itype][jtype] *
      (1.0/(denlj*denlj) - 1.0/denlj);
    if (rsq > cut_lj_innersq) {
      switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
        (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) / denom_lj;
      philj *= switch1;
    }
    eng += factor_lj*philj;
  }

  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairLJCharmmCoulLongSoft::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"lj14_1") == 0) return (void *) lj14_1;
  if (strcmp(str,"lj14_2") == 0) return (void *) lj14_2;
  if (strcmp(str,"lj14_3") == 0) return (void *) lj14_3;
  if (strcmp(str,"lj14_4") == 0) return (void *) lj14_4;

  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"sigma") == 0) return (void *) sigma;
  if (strcmp(str,"lambda") == 0) return (void *) lambda;

  dim = 0;
  if (strcmp(str,"implicit") == 0) return (void *) &implicit;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;

  return NULL;
}
