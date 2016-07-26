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
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_charmm_coul_long_14.h"
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

PairLJCharmmCoulLong14::PairLJCharmmCoulLong14(LAMMPS *lmp) : PairLJCharmmCoulLong(lmp)
{
  respa_enable = 1;
  ewaldflag = pppmflag = 1;
  ftable = NULL;
  implicit = 0;
  mix_flag = ARITHMETIC;
  writedata = 1;
  //allocated = 0;
}

/* ---------------------------------------------------------------------- */

PairLJCharmmCoulLong14::~PairLJCharmmCoulLong14()
{
  memory->destroy(elescale);
}

/* ---------------------------------------------------------------------- */

void PairLJCharmmCoulLong14::compute(int eflag, int vflag)
{
  int i,j,ii,jj,jo,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double **this_lj1,**this_lj2,**this_lj3,**this_lj4;
  double grij,expm2,prefactor,t,erfc;
  double philj,switch1,switch2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

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
      jo = j = jlist[jj];
      j &= NEIGHMASK; 
      jtype = type[j];
      factor_lj = special_lj[sbmask(jo)];
      factor_coul = special_coul[sbmask(jo)];
      if (sbmask(jo) == 3)  {
        this_lj1 = lj14_1;
        this_lj2 = lj14_2;
        this_lj3 = lj14_3;
        this_lj4 = lj14_4 ;        
        // if this is not true, then the 1-4 interaction are already pre-scaled
        if (!((lj3[itype][jtype]==lj14_3[itype][jtype])&&(lj4[itype][jtype]==lj14_4[itype][jtype])))
          factor_lj = 1.0;
      } else {
        this_lj1 = lj1;
        this_lj2 = lj2;
        this_lj3 = lj3;
        this_lj4 = lj4;
      }
      qtmp *= elescale[itype][jtype];     

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_bothsq) {
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
            prefactor = qqrd2e * qtmp*q[j]/r;
            forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
            if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
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

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;          
          forcelj = r6inv * (this_lj1[itype][jtype]*r6inv - this_lj2[itype][jtype]);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
              (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
            switch2 = 12.0*rsq * (cut_ljsq-rsq) *
              (rsq-cut_lj_innersq) * denom_lj_inv;
            philj = r6inv * (this_lj3[itype][jtype]*r6inv - this_lj4[itype][jtype]);
            forcelj = forcelj*switch1 + philj*switch2;
          }
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

          if (rsq < cut_ljsq) {
            evdwl = r6inv*(this_lj3[itype][jtype]*r6inv-this_lj4[itype][jtype]);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
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

void PairLJCharmmCoulLong14::compute_inner()
{
  int i,j,ii,jj,jo,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double **this_lj1,**this_lj2,**this_lj3,**this_lj4;
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
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      jo = j = jlist[jj];
      j &= NEIGHMASK; 
      jtype = type[j];
      factor_lj = special_lj[sbmask(jo)];
      factor_coul = special_coul[sbmask(jo)];
      if (sbmask(j) == 3)  {
        this_lj1 = lj14_1;
        this_lj2 = lj14_2;
        this_lj3 = lj14_3;
        this_lj4 = lj14_4 ;        
        // if this is not true, then the 1-4 interaction are already pre-scaled
        if (!((lj3[itype][jtype]==lj14_3[itype][jtype])&&(lj4[itype][jtype]==lj14_4[itype][jtype])))
          factor_lj = 1.0;
      } else {
        this_lj1 = lj1;
        this_lj2 = lj2;
        this_lj3 = lj3;
        this_lj4 = lj4;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_out_off_sq) {
        r2inv = 1.0/rsq;
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        r6inv = r2inv*r2inv*r2inv;
        jtype = type[j];
        forcelj = r6inv * (this_lj1[itype][jtype]*r6inv - this_lj2[itype][jtype]);

        fpair = (forcecoul + factor_lj*forcelj) * r2inv;

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

void PairLJCharmmCoulLong14::compute_middle()
{
  int i,j,ii,jj,jo,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,fpair;
  double rsq,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double **this_lj1,**this_lj2,**this_lj3,**this_lj4;
  double philj,switch1,switch2;
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
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      jo = j = jlist[jj];
      j &= NEIGHMASK; 
      jtype = type[j];
      factor_lj = special_lj[sbmask(jo)];
      factor_coul = special_coul[sbmask(jo)];
      if (sbmask(j) == 3)  {
        this_lj1 = lj14_1;
        this_lj2 = lj14_2;
        this_lj3 = lj14_3;
        this_lj4 = lj14_4 ;        
        // if this is not true, then the 1-4 interaction are already pre-scaled
        if (!((lj3[itype][jtype]==lj14_3[itype][jtype])&&(lj4[itype][jtype]==lj14_4[itype][jtype])))
          factor_lj = 1.0;
      } else {
        this_lj1 = lj1;
        this_lj2 = lj2;
        this_lj3 = lj3;
        this_lj4 = lj4;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_out_off_sq && rsq > cut_in_off_sq) {
        r2inv = 1.0/rsq;
        forcecoul = qqrd2e * qtmp*q[j]*sqrt(r2inv);
        if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*forcecoul;

        r6inv = r2inv*r2inv*r2inv;
        forcelj = r6inv * (this_lj1[itype][jtype]*r6inv - this_lj2[itype][jtype]);
        if (rsq > cut_lj_innersq) {
          switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
            (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
          switch2 = 12.0*rsq * (cut_ljsq-rsq) *
            (rsq-cut_lj_innersq) * denom_lj_inv;
          philj = r6inv * (this_lj3[itype][jtype]*r6inv - this_lj4[itype][jtype]);
          forcelj = forcelj*switch1 + philj*switch2;
        }

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

void PairLJCharmmCoulLong14::compute_outer(int eflag, int vflag)
{
  int i,j,ii,jj,jo,inum,jnum,itype,jtype,itable;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double r,r2inv,r6inv,forcecoul,forcelj,factor_coul,factor_lj;
  double **this_lj1,**this_lj2,**this_lj3,**this_lj4;
  double grij,expm2,prefactor,t,erfc;
  double philj,switch1,switch2;
  double rsw;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

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
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      jo = j = jlist[jj];
      j &= NEIGHMASK; 
      jtype = type[j];
      factor_lj = special_lj[sbmask(jo)];
      factor_coul = special_coul[sbmask(jo)];
      if (sbmask(j) == 3)  {
        this_lj1 = lj14_1;
        this_lj2 = lj14_2;
        this_lj3 = lj14_3;
        this_lj4 = lj14_4 ;                
        // if this is not true, then the 1-4 interaction are already pre-scaled
        if (!((lj3[itype][jtype]==lj14_3[itype][jtype])&&(lj4[itype][jtype]==lj14_4[itype][jtype])))
          factor_lj = 1.0;
      } else {
        this_lj1 = lj1;
        this_lj2 = lj2;
        this_lj3 = lj3;
        this_lj4 = lj4;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq < cut_bothsq) {
        r2inv = 1.0/rsq;

        if (rsq < cut_coulsq) {
          if (!ncoultablebits || rsq <= tabinnersq) {
            r = sqrt(rsq);
            grij = g_ewald * r;
            expm2 = exp(-grij*grij);
            t = 1.0 / (1.0 + EWALD_P*grij);
            erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
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

        if (rsq < cut_ljsq && rsq > cut_in_off_sq) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv * (this_lj1[itype][jtype]*r6inv - this_lj2[itype][jtype]);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
              (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
            switch2 = 12.0*rsq * (cut_ljsq-rsq) *
              (rsq-cut_lj_innersq) * denom_lj_inv;
            philj = r6inv * (this_lj3[itype][jtype]*r6inv - this_lj4[itype][jtype]);
            forcelj = forcelj*switch1 + philj*switch2;
          }
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

          if (rsq < cut_ljsq) {
            r6inv = r2inv*r2inv*r2inv;
            evdwl = r6inv*(this_lj3[itype][jtype]*r6inv-this_lj4[itype][jtype]);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
              evdwl *= switch1;
            }
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
            r6inv = r2inv*r2inv*r2inv;
            forcelj = r6inv * (this_lj1[itype][jtype]*r6inv - this_lj2[itype][jtype]);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
              switch2 = 12.0*rsq * (cut_ljsq-rsq) *
                (rsq-cut_lj_innersq) * denom_lj_inv;
              philj = r6inv * (this_lj3[itype][jtype]*r6inv - this_lj4[itype][jtype]);
              forcelj = forcelj*switch1 + philj*switch2;
            }
          } else if (rsq <= cut_in_on_sq) {
            forcelj = r6inv * (this_lj1[itype][jtype]*r6inv - this_lj2[itype][jtype]);
            if (rsq > cut_lj_innersq) {
              switch1 = (cut_ljsq-rsq) * (cut_ljsq-rsq) *
                (cut_ljsq + 2.0*rsq - 3.0*cut_lj_innersq) * denom_lj_inv;
              switch2 = 12.0*rsq * (cut_ljsq-rsq) *
                (rsq-cut_lj_innersq) * denom_lj_inv;
              philj = r6inv * (this_lj3[itype][jtype]*r6inv - this_lj4[itype][jtype]);
              forcelj = forcelj*switch1 + philj*switch2;
            }
          }

          fpair = (forcecoul + factor_lj*forcelj) * r2inv;
        }

        if (evflag) ev_tally(i,j,nlocal,newton_pair,
                             evdwl,ecoul,fpair,delx,dely,delz);
      }
    }
  }
}


/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCharmmCoulLong14::init_style()
{
  PairLJCharmmCoulLong::init_style();
  int n = atom->ntypes;

  memory->create(elescale,n+1,n+1,"pair:elescale");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++) {
      elescale[i][j] = 1.0;
      elescale[j][i] = 1.0;
    }
}

/* ---------------------------------------------------------------------- */

void *PairLJCharmmCoulLong14::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str,"lj14_1") == 0) return (void *) lj14_1;
  if (strcmp(str,"lj14_2") == 0) return (void *) lj14_2;
  if (strcmp(str,"lj14_3") == 0) return (void *) lj14_3;
  if (strcmp(str,"lj14_4") == 0) return (void *) lj14_4;
  if (strcmp(str,"epsilon") == 0) return (void *) epsilon;
  if (strcmp(str,"elescale") == 0) return (void *) elescale;

  dim = 0;
  if (strcmp(str,"implicit") == 0) return (void *) &implicit;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;

  return NULL;
}
