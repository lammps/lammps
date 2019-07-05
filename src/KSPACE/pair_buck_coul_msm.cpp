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

#include "pair_buck_coul_msm.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace MathConst;


/* ---------------------------------------------------------------------- */

PairBuckCoulMSM::PairBuckCoulMSM(LAMMPS *lmp) : PairBuckCoulLong(lmp)
{
  ewaldflag = pppmflag = 0;
  msmflag = 1;
  nmax = 0;
  ftmp = NULL;
}

/* ---------------------------------------------------------------------- */

PairBuckCoulMSM::~PairBuckCoulMSM()
{
  if (ftmp) memory->destroy(ftmp);
}

/* ---------------------------------------------------------------------- */

void PairBuckCoulMSM::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul,fpair,fcoul;
  double rsq,r2inv,r6inv,forcecoul,forcebuck,factor_coul,factor_lj;
  double egamma,fgamma,prefactor;
  double r,rexp;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int eflag_old = eflag;

  if (force->kspace->scalar_pressure_flag && vflag) {
    if (vflag > 2)
      error->all(FLERR,"Must use 'kspace_modify pressure/scalar no' "
        "to obtain per-atom virial with kspace_style MSM");

    if (atom->nmax > nmax) {
      if (ftmp) memory->destroy(ftmp);
      nmax = atom->nmax;
      memory->create(ftmp,nmax,3,"pair:ftmp");
    }
    memset(&ftmp[0][0],0,nmax*3*sizeof(double));

    // must switch on global energy computation if not already on

    if (eflag == 0 || eflag == 2) {
      eflag++;
    }
  }

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
        r2inv = 1.0/rsq;
        r = sqrt(rsq);

        if (rsq < cut_coulsq) {
          prefactor = qqrd2e * qtmp*q[j]/r;
          egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
          fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
          forcecoul = prefactor * fgamma;
          if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
        } else forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          r6inv = r2inv*r2inv*r2inv;
          rexp = exp(-r*rhoinv[itype][jtype]);
          forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
        } else forcebuck = 0.0;

        if (!(force->kspace->scalar_pressure_flag && vflag)) {
          fpair = (forcecoul + factor_lj*forcebuck) * r2inv;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }
        } else {
          // separate Buck and Coulombic forces

          fpair = (factor_lj*forcebuck) * r2inv;

          f[i][0] += delx*fpair;
          f[i][1] += dely*fpair;
          f[i][2] += delz*fpair;
          if (newton_pair || j < nlocal) {
            f[j][0] -= delx*fpair;
            f[j][1] -= dely*fpair;
            f[j][2] -= delz*fpair;
          }

          fcoul = (forcecoul) * r2inv;

          ftmp[i][0] += delx*fcoul;
          ftmp[i][1] += dely*fcoul;
          ftmp[i][2] += delz*fcoul;
          if (newton_pair || j < nlocal) {
            ftmp[j][0] -= delx*fcoul;
            ftmp[j][1] -= dely*fcoul;
            ftmp[j][2] -= delz*fcoul;
          }
        }

        if (eflag) {
          if (rsq < cut_coulsq) {
            ecoul = prefactor*egamma;
            if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
          } else ecoul = 0.0;
          if (eflag_old && rsq < cut_ljsq[itype][jtype]) {
            evdwl = a[itype][jtype]*rexp - c[itype][jtype]*r6inv -
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

  if (force->kspace->scalar_pressure_flag && vflag) {
    for (i = 0; i < 3; i++) virial[i] += force->pair->eng_coul/3.0;
    for (int i = 0; i < nmax; i++) {
      f[i][0] += ftmp[i][0];
      f[i][1] += ftmp[i][1];
      f[i][2] += ftmp[i][2];
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairBuckCoulMSM::single(int i, int j, int itype, int jtype,
                                double rsq,
                                double factor_coul, double factor_lj,
                                double &fforce)
{
  double r2inv,r6inv,r,rexp,egamma,fgamma,prefactor;
  double forcecoul,forcebuck,phicoul,phibuck;

  r2inv = 1.0/rsq;
  if (rsq < cut_coulsq) {
    r = sqrt(rsq);
    prefactor = force->qqrd2e * atom->q[i]*atom->q[j]/r;
    egamma = 1.0 - (r/cut_coul)*force->kspace->gamma(r/cut_coul);
    fgamma = 1.0 + (rsq/cut_coulsq)*force->kspace->dgamma(r/cut_coul);
    forcecoul = prefactor * fgamma;
    if (factor_coul < 1.0) forcecoul -= (1.0-factor_coul)*prefactor;
  } else forcecoul = 0.0;
  if (rsq < cut_ljsq[itype][jtype]) {
    r6inv = r2inv*r2inv*r2inv;
    r = sqrt(rsq);
    rexp = exp(-r*rhoinv[itype][jtype]);
    forcebuck = buck1[itype][jtype]*r*rexp - buck2[itype][jtype]*r6inv;
  } else forcebuck = 0.0;
  fforce = (forcecoul + factor_lj*forcebuck) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq) {
    phicoul = prefactor*egamma;
    if (factor_coul < 1.0) phicoul -= (1.0-factor_coul)*prefactor;
    eng += phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) {
    phibuck = a[itype][jtype]*rexp - c[itype][jtype]*r6inv -
      offset[itype][jtype];
    eng += factor_lj*phibuck;
  }
  return eng;
}

/* ---------------------------------------------------------------------- */

void *PairBuckCoulMSM::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  return NULL;
}
