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
   Contributing author: Hendrik Heenen (hendrik.heenen@mytum.de)
------------------------------------------------------------------------- */

#include "pair_born_coul_long_cs.h"
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

PairBornCoulLongCS::PairBornCoulLongCS(LAMMPS *lmp) : PairBornCoulLong(lmp)
{
  ewaldflag = pppmflag = 1;
  single_enable = 0; // TODO: single function does not match compute
  ftable = nullptr;
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

void PairBornCoulLongCS::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itable,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair;
  double fraction,table;
  double rsq,r2inv,r6inv,forcecoul,forceborn,factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc,u;
  double r,rexp;
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
              // used approximation functions valid
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
          r = sqrt(rsq);
          r6inv = r2inv*r2inv*r2inv;
          rexp = exp((sigma[itype][jtype]-r)*rhoinv[itype][jtype]);
          forceborn = born1[itype][jtype]*r*rexp - born2[itype][jtype]*r6inv
            + born3[itype][jtype]*r2inv*r6inv;
        } else forceborn = 0.0;

        fpair = (forcecoul + factor_lj*forceborn) * r2inv;

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
            evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv
              + d[itype][jtype]*r6inv*r2inv - offset[itype][jtype];
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
