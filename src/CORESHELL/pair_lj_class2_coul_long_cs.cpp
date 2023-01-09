// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_class2_coul_long_cs.h"
#include <cmath>
#include "atom.h"
#include "force.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

#define EWALD_F   1.12837917
#define EWALD_P   9.95473818e-1
#define B0       -0.1335096380159268
#define B1       -2.57839507e-1
#define B2       -1.37203639e-1
#define B3       -8.88822059e-3
#define B4       -5.80844129e-3
#define B5        1.14652755e-1

#define EPSILON 1.0e-20
#define EPS_EWALD 1.0e-6
#define EPS_EWALD_SQR 1.0e-12

/* ---------------------------------------------------------------------- */

PairLJClass2CoulLongCS::PairLJClass2CoulLongCS(LAMMPS *lmp) : PairLJClass2CoulLong(lmp)
{
  ewaldflag = pppmflag = 1;
  respa_enable = 0;  // TODO: r-RESPA handling is inconsistent and thus disabled until fixed
  single_enable = 0; // TODO: single function does not match compute
  writedata = 1;
  ftable = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJClass2CoulLongCS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double rsq,r,rinv,r2inv,r3inv,r6inv,forcecoul,forcelj;
  double grij,expm2,prefactor,t,erfc,u;
  double factor_coul,factor_lj;
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

      if (rsq < cutsq[itype][jtype]) {
        rsq += EPSILON; // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
        r2inv = 1.0/rsq;
        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            prefactor = qqrd2e * qtmp*q[j];
            if (factor_coul < 1.0) {
              // When bonded parts are being calculated a minimal distance (EPS_EWALD)
              // has to be added to the prefactor and erfc in order to make the
              // used approximation functions for the Ewald correction valid
              grij = g_ewald * (r+EPS_EWALD);
              expm2 = exp(-grij*grij);
              t = 1.0 / (1.0 + EWALD_P*grij);
              u = 1.0 - t;
              erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
              prefactor /= (r+EPS_EWALD);
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - (1.0-factor_coul));
              // Additionally r2inv needs to be accordingly modified since the later
              // scaling of the overall force shall be consistent
              r2inv = 1.0/(rsq + EPS_EWALD_SQR);
            } else {
              grij = g_ewald * r;
              expm2 = exp(-grij*grij);
              t = 1.0 / (1.0 + EWALD_P*grij);
              u = 1.0 - t;
              erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
              prefactor /= r;
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          rinv = sqrt(r2inv);
          r3inv = r2inv*rinv;
          r6inv = r3inv*r3inv;
          forcelj = r6inv * (lj1[itype][jtype]*r3inv - lj2[itype][jtype]);
        } else forcelj = 0.0;

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;

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
            if (!ncoultablebits || rsq <= tabinnersq)
              ecoul = prefactor*erfc;
            else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
            }
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            evdwl = r6inv*(lj3[itype][jtype]*r3inv-lj4[itype][jtype]) -
              offset[itype][jtype];
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

void PairLJClass2CoulLongCS::compute_inner()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,rinv,r2inv,r3inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double rsw;
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
        rsq += EPSILON; // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
        r2inv = 1.0/rsq;
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        jtype = type[j];
        if (rsq < cut_ljsq[itype][jtype]) {
          rinv = sqrt(r2inv);
          r3inv = r2inv*rinv;
          r6inv = r3inv*r3inv;
          forcelj = r6inv * (lj1[itype][jtype]*r3inv - lj2[itype][jtype]);
        } else forcelj = 0.0;

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;
        if (rsq > cut_out_on_sq) {
          rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
          fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
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

void PairLJClass2CoulLongCS::compute_middle()
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,rinv,r2inv,r3inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double rsw;
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
        r2inv = 1.0/rsq;
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        jtype = type[j];
        if (rsq < cut_ljsq[itype][jtype]) {
          rinv = sqrt(r2inv);
          r3inv = r2inv*rinv;
          r6inv = r3inv*r3inv;
          forcelj = r6inv * (lj1[itype][jtype]*r3inv - lj2[itype][jtype]);
        } else forcelj = 0.0;

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;
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

void PairLJClass2CoulLongCS::compute_outer(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,rinv,r2inv,r3inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc,u;
  double rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

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
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            u = 1. - t;
            erfc = t * (1.+u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2 - 1.0);
            if (rsq > cut_in_off_sq) {
              if (rsq < cut_in_on_sq) {
                rsw = (r - cut_in_off)/cut_in_diff;
                forcecoul += prefactor*rsw*rsw*(3.0 - 2.0*rsw);
                if (factor_coul < 1.0)
                  forcecoul -=
                    (1.0-factor_coul)*prefactor*rsw*rsw*(3.0 - 2.0*rsw);
              } else {
                forcecoul += prefactor;
                if (factor_coul < 1.0)
                  forcecoul -= (1.0-factor_coul)*prefactor;
              }
            }
          } else {
            union_int_float_t rsq_lookup;
            rsq_lookup.f = rsq;
            itable = rsq_lookup.i & ncoulmask;
            itable >>= ncoulshiftbits;
            fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
            table = ftable[itable] + fraction*dftable[itable];
            forcecoul = qtmp*q[j] * table;
            if (factor_coul < 1.0) {
              table = ctable[itable] + fraction*dctable[itable];
              prefactor = qtmp*q[j] * table;
              forcecoul -= (1.0-factor_coul)*prefactor;
            }
          }
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype] && rsq > cut_in_off_sq) {
          rinv = sqrt(r2inv);
          r3inv = r2inv*rinv;
          r6inv = r3inv*r3inv;
          forcelj = r6inv * (lj1[itype][jtype]*r3inv - lj2[itype][jtype]);
          if (rsq < cut_in_on_sq) {
            rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
            forcelj *= rsw*rsw*(3.0 - 2.0*rsw);
          }
        } else forcelj = 0.0;

        fpair = (forcecoul + forcelj) * r2inv;

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
            if (!ncoultablebits || rsq <= tabinnersq) {
              ecoul = prefactor*erfc;
              if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
            } else {
              table = etable[itable] + fraction*detable[itable];
              ecoul = qtmp*q[j] * table;
              if (factor_coul < 1.0) {
                table = ptable[itable] + fraction*dptable[itable];
                prefactor = qtmp*q[j] * table;
                ecoul -= (1.0-factor_coul)*prefactor;
              }
            }
          } else ecoul = 0.0;

          if (rsq < cut_ljsq[itype][jtype]) {
            rinv = sqrt(r2inv);
            r3inv = r2inv*rinv;
            r6inv = r3inv*r3inv;
            evdwl = r6inv*(lj3[itype][jtype]*r3inv-lj4[itype][jtype]) -
              offset[itype][jtype];
            evdwl *= factor_lj;
          } else evdwl = 0.0;
        }

        if (vflag) {
          if (rsq < cut_coulsq) {
            if (!ncoultablebits || rsq <= tabinnersq) {
              forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
              if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
            } else {
              table = vtable[itable] + fraction*dvtable[itable];
              forcecoul = qtmp*q[j] * table;
              if (factor_coul < 1.0) {
                table = ptable[itable] + fraction*dptable[itable];
                prefactor = qtmp*q[j] * table;
                forcecoul -= (1.0-factor_coul)*prefactor;
              }
            }
          } else forcecoul = 0.0;

          if (rsq <= cut_in_off_sq) {
            rinv = sqrt(r2inv);
            r3inv = r2inv*rinv;
            r6inv = r3inv*r3inv;
            forcelj = r6inv * (lj1[itype][jtype]*r3inv - lj2[itype][jtype]);
          } else if (rsq <= cut_in_on_sq) {
            rinv = sqrt(r2inv);
            r3inv = r2inv*rinv;
            r6inv = r3inv*r3inv;
            forcelj = r6inv * (lj1[itype][jtype]*r3inv - lj2[itype][jtype]);
          }
          fpair = (forcecoul + factor_lj*forcelj) * r2inv;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}

